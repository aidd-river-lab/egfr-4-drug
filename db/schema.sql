-- EGFR C797S 四代抑制剂研发管线 —— 数据库schema
--
-- 设计原则(和 core/ 下所有模块一致)：
--   1. 字段名尽量直接照抄各 core/*.py 里对应 dataclass 的字段名，方便 ORM/脚本直接映射，
--      不要在这里发明一套新命名。
--   2. 允许大量字段为 NULL —— 很多列(FEP、MM-GBSA、共价对接、湿实验ADME)在本环境里
--      还没有真实计算/测量，表结构要如实反映"有schema但没数据"，不能为了让表"看起来满"
--      而给默认值。
--   3. 用 SQLite 方言写(本地单机就能跑，sqlite3 schema.sql 直接建库)；如果以后换
--      Postgres，主要改动是 AUTOINCREMENT->SERIAL/IDENTITY，CHECK约束基本通用。
--
-- 建库: sqlite3 egfr4.db < schema.sql

PRAGMA foreign_keys = ON;

-- ============================================================
-- 环节1：标准化后的化合物主表(core/standardize.py StandardizeResult.curated)
-- ============================================================
CREATE TABLE compounds (
    compound_id         TEXT PRIMARY KEY,           -- 内部编号，比如 DEMO-001
    smiles_canonical    TEXT NOT NULL,
    inchikey            TEXT NOT NULL,
    inchikey_skeleton   TEXT NOT NULL,               -- 前14位，用于忽略立体/电荷的去重查重
    source              TEXT,                        -- 来源：enumeration / literature / vendor / wet_lab
    scaffold_id         TEXT,                        -- 关联 scaffolds.yaml 里的scaffold名字
    created_at          TEXT NOT NULL DEFAULT (datetime('now')),
    UNIQUE (inchikey)
);

CREATE INDEX idx_compounds_inchikey_skeleton ON compounds (inchikey_skeleton);

-- 标准化时被隔离的记录(core/standardize.py StandardizeResult.quarantine)，
-- 不进 compounds 表，但要留痕方便回查原始数据有什么问题
CREATE TABLE standardize_quarantine (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    raw_smiles          TEXT NOT NULL,
    reason              TEXT NOT NULL,
    quarantined_at      TEXT NOT NULL DEFAULT (datetime('now'))
);

-- ============================================================
-- 环节2：结构系综(core/structures.py StructureEnsembleMember)
-- ============================================================
CREATE TABLE structure_ensemble (
    ensemble_id             TEXT PRIMARY KEY,          -- 比如 del19_C797S_cluster1
    genotype                TEXT NOT NULL,             -- 对应 target_profile.yaml 的基因型id
    source_pdb              TEXT NOT NULL,
    mutation_method         TEXT,                      -- Rosetta residue mutation / none(实验结构)
    md_trajectory_id        TEXT,
    cluster_occupancy_pct   REAL,
    prep_date               TEXT,
    forcefield_version      TEXT,
    local_path              TEXT,
    notes                   TEXT DEFAULT ''
);

CREATE INDEX idx_structure_ensemble_genotype ON structure_ensemble (genotype);

-- ============================================================
-- 环节3：枚举出来的设计分子，和具体scaffold/R基团组合关联
-- ============================================================
CREATE TABLE designed_molecules (
    compound_id         TEXT PRIMARY KEY REFERENCES compounds(compound_id),
    scaffold_id         TEXT NOT NULL,
    r_group_assignment  TEXT NOT NULL,   -- JSON，比如 {"R1":"...","R2":"...","R3":"..."}，对应pocket_regions.yaml的fragment id
    enumeration_batch   TEXT,
    created_at          TEXT NOT NULL DEFAULT (datetime('now'))
);

-- ============================================================
-- 环节4：多层打分漏斗 L0-L4(core/docking.py DockingResult, core/mmgbsa.py
-- MMGBSAResult, core/fep.py FEPPerturbation/FEPMapResult 的落库形式)
-- 一个化合物在一个结构系综成员上，每一层最多一条记录
-- ============================================================
CREATE TABLE funnel_scores (
    id                      INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id             TEXT NOT NULL REFERENCES compounds(compound_id),
    ensemble_id             TEXT NOT NULL REFERENCES structure_ensemble(ensemble_id),
    funnel_level            TEXT NOT NULL CHECK (funnel_level IN ('L0', 'L1', 'L2', 'L3', 'L4')),
    method                  TEXT NOT NULL,             -- pharmacophore / vina_docking / refined_docking / md_mmgbsa / fep
    score_value             REAL,                      -- 统一单位 kcal/mol(pharmacophore层除外，记匹配分)
    error_estimate_kcal_mol REAL,                      -- 对应 pipeline.yaml funnel.<L>.error_kcal_mol
    n_poses_or_frames       INTEGER,
    passed_filter           INTEGER CHECK (passed_filter IN (0, 1)),  -- 是否通过该层pass_criteria
    raw_output_path         TEXT,                      -- pdbqt/轨迹/原始输出文件路径
    computed_at             TEXT NOT NULL DEFAULT (datetime('now')),
    UNIQUE (compound_id, ensemble_id, funnel_level)
);

CREATE INDEX idx_funnel_scores_compound ON funnel_scores (compound_id);

-- L3 的位姿/氢键质检明细(pipeline.yaml funnel.L3.pass_criteria: rmsd/hinge_hbond/target_anchor)
CREATE TABLE md_stability_qc (
    funnel_score_id         INTEGER PRIMARY KEY REFERENCES funnel_scores(id),
    rmsd_angstrom           REAL,
    hinge_hbond_occupancy_pct REAL,
    target_anchor_occupancy_pct REAL,
    qc_pass                 INTEGER CHECK (qc_pass IN (0, 1))
);

-- L4 FEP map的质控(core/fep.py FEPMapResult.qc_pass逻辑的落库)
CREATE TABLE fep_maps (
    fep_map_id          TEXT PRIMARY KEY,
    genotype            TEXT NOT NULL,
    mue_vs_experiment   REAL,            -- core/fep.py mue_vs_experiment()
    n_anchor_pairs      INTEGER NOT NULL DEFAULT 0,
    qc_pass             INTEGER CHECK (qc_pass IN (0, 1)),
    rejected_reason      TEXT,           -- qc_pass=0时必须填，比如"MUE超过1.5kcal/mol阈值"
    computed_at         TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE fep_perturbations (
    id                          INTEGER PRIMARY KEY AUTOINCREMENT,
    fep_map_id                  TEXT NOT NULL REFERENCES fep_maps(fep_map_id),
    ligand_a                    TEXT NOT NULL REFERENCES compounds(compound_id),
    ligand_b                    TEXT NOT NULL REFERENCES compounds(compound_id),
    predicted_ddg_kcal_mol      REAL NOT NULL,
    uncertainty_kcal_mol        REAL NOT NULL,
    experimental_ddg_kcal_mol   REAL         -- 非NULL表示这对是实验锚点
);

-- ============================================================
-- 环节4.3：共价对接(core/covalent.py CovalentDockingResult)
-- ============================================================
CREATE TABLE covalent_docking (
    compound_id                                TEXT NOT NULL REFERENCES compounds(compound_id),
    ensemble_id                                TEXT NOT NULL REFERENCES structure_ensemble(ensemble_id),
    nucleophile_residue                        TEXT NOT NULL,
    warhead_carbon_to_nucleophile_distance_a   REAL,
    attack_angle_degree                        REAL,
    attack_geometry_ok                         INTEGER CHECK (attack_geometry_ok IN (0, 1)),  -- Bürgi-Dunitz判据
    ki_nm                                      REAL,          -- 可逆亲和力，决定选择性
    kinact_per_second                          REAL,          -- 成键速率
    gsh_half_life_hours                        REAL,          -- 弹头本征反应性
    PRIMARY KEY (compound_id, ensemble_id)
    -- kinact_over_ki 不落库，是只读计算属性(kinact/KI)，查询时现算，避免和源字段不同步
);

-- ============================================================
-- 环节5：选择性引擎(core/selectivity.py SelectivityResult)
-- ============================================================
CREATE TABLE selectivity_results (
    id                      INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id             TEXT NOT NULL REFERENCES compounds(compound_id),
    genotype                TEXT NOT NULL,           -- 目标突变体，比如 del19_C797S
    ddg_target_kcal_mol     REAL,
    ddg_wt_kcal_mol         REAL,
    dddg_kcal_mol           REAL,
    fold_selectivity_pred   REAL,
    uncertainty_kcal_mol    REAL,
    method                  TEXT NOT NULL CHECK (method IN ('fep', 'mmgbsa', 'ml', 'docking')),
    structures_used         TEXT,                    -- JSON数组，ensemble_id列表
    computed_at             TEXT NOT NULL DEFAULT (datetime('now')),
    UNIQUE (compound_id, genotype, method)
);

-- ============================================================
-- 环节6：ADMET描述符 + 结构警示(core/admet.py)
-- ============================================================
CREATE TABLE admet_descriptors (
    compound_id             TEXT PRIMARY KEY REFERENCES compounds(compound_id),
    mw                      REAL,
    clogp                   REAL,
    clogd_approx            REAL,          -- 近似值，见 core/admet.py 顶部声明
    tpsa                    REAL,
    hbd                     INTEGER,
    hba                     INTEGER,
    fsp3                    REAL,
    has_basic_aliphatic_amine INTEGER CHECK (has_basic_aliphatic_amine IN (0, 1)),
    pka_proxy               TEXT,          -- 子结构代理标记，不是数值pKa
    cns_mpo                 REAL,
    cns_mpo_pass_4_0        INTEGER CHECK (cns_mpo_pass_4_0 IN (0, 1)),
    herg_risk_proxy         TEXT CHECK (herg_risk_proxy IN ('low', 'medium', 'high')),
    computed_at             TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE structural_alert_hits (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id     TEXT NOT NULL REFERENCES compounds(compound_id),
    rule_name       TEXT NOT NULL,          -- 对应 structural_alerts.yaml 里的key
    severity        TEXT NOT NULL,
    action          TEXT NOT NULL CHECK (action IN ('PASS', 'OPTIMIZE_NEEDED', 'REJECT')),
    advice          TEXT,
    hit_at          TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_structural_alert_hits_compound ON structural_alert_hits (compound_id);

-- ============================================================
-- 环节7：合成可及性(core/synthesis.py SynthesisRecord)
-- ============================================================
CREATE TABLE synthesis_records (
    compound_id                 TEXT PRIMARY KEY REFERENCES compounds(compound_id),
    sa_score                    REAL,
    n_steps                     INTEGER,
    longest_linear_sequence     INTEGER,
    key_bb_cas                  TEXT,
    bb_lead_time_days           INTEGER,
    estimated_cost_per_20mg     REAL,
    chiral_resolution_needed    INTEGER CHECK (chiral_resolution_needed IN (0, 1)),
    scale_up_risk               TEXT,
    fto_status                  TEXT        -- in_claim / near_claim / clear / unknown_no_real_patent_data
);

-- ============================================================
-- 环节8：多目标决策(core/decide.py 的输出留痕，不存算法本身)
-- ============================================================
CREATE TABLE decision_rounds (
    round_id            TEXT PRIMARY KEY,         -- 比如 "2026-W40"
    round_date          TEXT NOT NULL,
    batch_size          INTEGER NOT NULL,
    quota_json          TEXT,                      -- JSON，记录当轮用的 exploit/explore/hypothesis/control 配额
    notes               TEXT DEFAULT ''
);

CREATE TABLE decision_batch_members (
    round_id            TEXT NOT NULL REFERENCES decision_rounds(round_id),
    compound_id         TEXT NOT NULL REFERENCES compounds(compound_id),
    batch_role          TEXT NOT NULL CHECK (batch_role IN ('exploit', 'explore', 'hypothesis_test', 'control')),
    is_pareto_optimal   INTEGER CHECK (is_pareto_optimal IN (0, 1)),
    desirability_score  REAL,
    hypothesis_label    TEXT,          -- hypothesis_test角色时，药化指定的具体假设描述(自由文本)
    PRIMARY KEY (round_id, compound_id)
);

-- ============================================================
-- 环节9：DMTA 回传数据(core/feedback.py WetLabResult)
-- ============================================================
CREATE TABLE wet_lab_results (
    id                      INTEGER PRIMARY KEY AUTOINCREMENT,
    -- 身份
    compound_id             TEXT NOT NULL REFERENCES compounds(compound_id),
    batch_id                TEXT NOT NULL,
    smiles_as_made          TEXT NOT NULL,     -- 实际做出来的结构，可能和设计的不同，不可省略
    purity_pct              REAL NOT NULL CHECK (purity_pct BETWEEN 0 AND 100),
    chirality_confirmed     INTEGER NOT NULL CHECK (chirality_confirmed IN (0, 1)),
    -- 活性
    genotype                TEXT NOT NULL,
    assay_format            TEXT NOT NULL CHECK (assay_format IN ('enzymatic', 'cellular', 'nanobret')),
    atp_conc_um              REAL,
    readout                 TEXT NOT NULL,
    value                   REAL NOT NULL,
    unit                    TEXT NOT NULL,
    censoring               TEXT NOT NULL DEFAULT 'point' CHECK (censoring IN ('point', 'less_than', 'greater_than')),
    n_replicates            INTEGER NOT NULL CHECK (n_replicates >= 1),
    log_sd                  REAL,
    assay_date              TEXT NOT NULL,
    plate_id                TEXT,
    -- 选择性
    wt_value                REAL,
    fold_selectivity        REAL,
    kinase_panel_file       TEXT,
    -- ADME
    sol_ph74                REAL,
    hlm_clint               REAL,
    rlm_clint               REAL,
    ppb_pct                 REAL,
    caco2_papp              REAL,
    mdr1_er                 REAL,
    herg_ic50_um            REAL,
    gsh_adduct_detected     INTEGER CHECK (gsh_adduct_detected IN (0, 1)),
    -- 体内
    species                 TEXT,
    dose_mg_per_kg          REAL,
    route                   TEXT,
    auc                     REAL,
    cmax                    REAL,
    t_half_hours            REAL,
    f_pct                   REAL,
    kp_uu                   REAL,
    -- 元数据
    operator                TEXT NOT NULL,
    protocol_version        TEXT NOT NULL,
    qc_status               TEXT NOT NULL DEFAULT 'pending' CHECK (qc_status IN ('pass', 'fail', 'pending')),
    comments                TEXT DEFAULT '',
    created_at              TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_wet_lab_results_compound ON wet_lab_results (compound_id);
CREATE INDEX idx_wet_lab_results_genotype ON wet_lab_results (genotype);

-- 假阳性归因(环节9.2第2条，标签必须结构化进数据库，不能自由文本——季度要统计)
CREATE TABLE false_positive_attributions (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id         TEXT NOT NULL REFERENCES compounds(compound_id),
    funnel_score_id     INTEGER REFERENCES funnel_scores(id),   -- 具体哪一层的预测被判定为假阳性
    attribution_label   TEXT NOT NULL CHECK (attribution_label IN (
        'pose_error', 'protonation_state_error', 'conformer_selection_error',
        'desolvation_underestimated', 'cell_permeability_issue',
        'metabolic_instability', 'compound_degradation_or_purity'
    )),
    attributed_by        TEXT NOT NULL,
    attributed_at        TEXT NOT NULL DEFAULT (datetime('now')),
    notes                TEXT DEFAULT ''
);

-- 规则治理：precision追踪(环节9.3，core/feedback.py rule_precision_tracking的落库，
-- 决定一条结构警示规则能不能从warn模式升级成reject)
CREATE TABLE rule_precision_history (
    id                          INTEGER PRIMARY KEY AUTOINCREMENT,
    rule_name                   TEXT NOT NULL,          -- 对应 structural_alerts.yaml 的key
    review_date                 TEXT NOT NULL,
    n_hits                      INTEGER NOT NULL,
    n_true_positive             INTEGER NOT NULL,
    precision                   REAL NOT NULL,
    n_independent_compounds     INTEGER NOT NULL,
    eligible_for_reject_mode    INTEGER NOT NULL CHECK (eligible_for_reject_mode IN (0, 1)),
    current_mode                TEXT NOT NULL CHECK (current_mode IN ('warn', 'reject'))
);

-- 前瞻验证时间序列(环节9.2第4条：模型重训后，每一轮的整体预测质量要能画成时间序列)
CREATE TABLE prospective_validation_rounds (
    round_id            TEXT PRIMARY KEY,
    round_date          TEXT NOT NULL,
    model_version       TEXT NOT NULL,
    n_compounds         INTEGER NOT NULL,
    spearman_rho        REAL,
    mue                 REAL,
    notes               TEXT DEFAULT ''
);
