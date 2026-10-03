-- EGFR C797S 第四代抑制剂研发管线 —— 数据库schema(MySQL 8.0+版本)
--
-- 这是 db/schema.sql(SQLite版) 的MySQL方言转换，表结构、字段名、字段含义完全一致，
-- 只改了类型系统和少数MySQL特有的写法。两份schema保持同步维护：SQLite版继续留着
-- 给不想启动MySQL server的场景做快速本地校验，MySQL版是生产环境的实际落地版本。
--
-- 2026-10 更新：compounds表加了route_id列，三条路线(routes/route_a_shp2_sos1、
-- routes/route_b_degrader、routes/route_c_4th_gen_tki)共用这一个数据库，唯一性
-- 约束也相应从UNIQUE(inchikey)改成UNIQUE(route_id, inchikey)，详见 db/schema.sql 顶部说明。
--
-- 转换时做的实质性调整(不是照抄改改类型名这么简单)：
--   1. 所有"身份类"文本字段(compound_id/ensemble_id/rule_name等)从SQLite的裸TEXT
--      改成有限长度的VARCHAR —— MySQL的InnoDB给TEXT/BLOB建索引必须指定前缀长度，
--      而这些字段几乎全部需要做主键/外键/UNIQUE/普通索引，裸TEXT在MySQL里用不了。
--      长度统一给得比观察到的最长真实值宽松一倍以上，省得以后遇到截断报错。
--   2. r_group_assignment / structures_used / quota_json 三个"JSON文本"字段
--      改成MySQL原生JSON类型(5.7.8+)，可以直接用 JSON_EXTRACT 等函数查，
--      比SQLite版本里"存JSON字符串但类型是TEXT"更好用。
--   3. CHECK约束需要 MySQL >= 8.0.16 才会真正生效(之前版本会被解析但不强制执行)，
--      本文件假设你的MySQL版本满足这个要求；达不到的话CHECK不报错但也不保护数据，
--      建议升级或者把关键约束前移到应用层。
--   4. 列名 `precision` 改成 `precision_value`——PRECISION是MySQL的保留字
--      (用在DOUBLE PRECISION这类类型声明里)，直接拿来做列名会出语法错误。
--      对应地，core/feedback.py::rule_precision_tracking() 返回的dict key仍然是
--      "precision"(Python端没有保留字冲突)，写入数据库时映射到这一列即可。
--   5. 时间戳统一用 DATETIME DEFAULT CURRENT_TIMESTAMP，纯日期字段(assay_date/
--      round_date/review_date/prep_date)用 DATE，不是全都用字符串占位。
--   6. 引擎统一 InnoDB(外键约束需要)，字符集统一 utf8mb4(存中文注释/建议文本，
--      utf8mb4而不是utf8是因为MySQL的"utf8"历史上是阉割版，不能存完整的emoji等
--      4字节字符，现在新建库没有理由不用utf8mb4)。
--
-- 建库：
--   mysql -u root -p < db/schema_mysql.sql
-- 这条命令会自动建一个叫 egfr4 的schema并切换进去，不需要提前手动建库。

CREATE DATABASE IF NOT EXISTS egfr4
    CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
USE egfr4;

-- ============================================================
-- 环节1：标准化后的化合物主表(core/standardize.py StandardizeResult.curated)
-- ============================================================
CREATE TABLE compounds (
    compound_id         VARCHAR(64) PRIMARY KEY,
    route_id            VARCHAR(32) NOT NULL DEFAULT 'route_c_4th_gen_tki',  -- route_a_shp2_sos1 / route_b_degrader / route_c_4th_gen_tki
    smiles_canonical    VARCHAR(2000) NOT NULL,
    inchikey            VARCHAR(32) NOT NULL,
    inchikey_skeleton   VARCHAR(20) NOT NULL,            -- 前14位，忽略立体/电荷的去重查重
    `source`              VARCHAR(32),                      -- enumeration / literature / vendor / wet_lab
    scaffold_id         VARCHAR(64),                       -- 关联该路线rgroup_libraries/scaffolds.yaml里的scaffold名字(配置文件,非DB外键)
    created_at          DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    UNIQUE KEY uq_compounds_route_inchikey (route_id, inchikey)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_compounds_inchikey_skeleton ON compounds (inchikey_skeleton);
CREATE INDEX idx_compounds_route_id ON compounds (route_id);

-- 标准化时被隔离的记录(core/standardize.py StandardizeResult.quarantine)
CREATE TABLE standardize_quarantine (
    id                  INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    raw_smiles          VARCHAR(2000) NOT NULL,
    reason              TEXT NOT NULL,
    quarantined_at      DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================
-- 环节2：结构系综(core/structures.py StructureEnsembleMember)
-- ============================================================
CREATE TABLE structure_ensemble (
    ensemble_id             VARCHAR(64) PRIMARY KEY,      -- 比如 del19_C797S_cluster1
    genotype                VARCHAR(64) NOT NULL,
    source_pdb              VARCHAR(16) NOT NULL,
    mutation_method         VARCHAR(255),                  -- "Rosetta residue mutation" / "none(实验结构)"
    md_trajectory_id        VARCHAR(64),
    cluster_occupancy_pct   DOUBLE,
    prep_date               DATE,
    forcefield_version      VARCHAR(64),
    local_path              VARCHAR(500),
    notes                   TEXT DEFAULT ''
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_structure_ensemble_genotype ON structure_ensemble (genotype);

-- ============================================================
-- 环节3：枚举出来的设计分子，和具体scaffold/R基团组合关联
-- ============================================================
CREATE TABLE designed_molecules (
    compound_id         VARCHAR(64) PRIMARY KEY,
    scaffold_id         VARCHAR(64) NOT NULL,
    r_group_assignment  JSON NOT NULL,                     -- {"R1":"...","R2":"...","R3":"..."}，对应pocket_regions.yaml的fragment id
    enumeration_batch   VARCHAR(64),
    created_at          DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_designed_molecules_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================
-- 环节4：多层打分漏斗 L0-L4
-- ============================================================
CREATE TABLE funnel_scores (
    id                      INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    compound_id             VARCHAR(64) NOT NULL,
    ensemble_id             VARCHAR(64) NOT NULL,
    funnel_level            VARCHAR(2) NOT NULL,
    `method`                  VARCHAR(64) NOT NULL,          -- pharmacophore / vina_docking / refined_docking / md_mmgbsa / fep
    score_value             DOUBLE,                         -- 统一单位kcal/mol(pharmacophore层除外，记匹配分)
    error_estimate_kcal_mol DOUBLE,                         -- 对应 pipeline.yaml funnel.<L>.error_kcal_mol
    n_poses_or_frames       INT,
    passed_filter           TINYINT(1),                     -- 是否通过该层pass_criteria
    raw_output_path         VARCHAR(500),
    computed_at             DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_funnel_scores_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT fk_funnel_scores_ensemble
        FOREIGN KEY (ensemble_id) REFERENCES structure_ensemble(ensemble_id),
    CONSTRAINT chk_funnel_scores_level CHECK (funnel_level IN ('L0', 'L1', 'L2', 'L3', 'L4')),
    CONSTRAINT chk_funnel_scores_passed CHECK (passed_filter IN (0, 1)),
    UNIQUE KEY uq_funnel_scores (compound_id, ensemble_id, funnel_level)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_funnel_scores_compound ON funnel_scores (compound_id);

-- L3 的位姿/氢键质检明细
CREATE TABLE md_stability_qc (
    funnel_score_id             INT UNSIGNED PRIMARY KEY,
    rmsd_angstrom                DOUBLE,
    hinge_hbond_occupancy_pct    DOUBLE,
    target_anchor_occupancy_pct  DOUBLE,
    qc_pass                      TINYINT(1),
    CONSTRAINT fk_md_stability_qc_funnel
        FOREIGN KEY (funnel_score_id) REFERENCES funnel_scores(id),
    CONSTRAINT chk_md_stability_qc_pass CHECK (qc_pass IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- L4 FEP map的质控(core/fep.py FEPMapResult.qc_pass逻辑的落库)
CREATE TABLE fep_maps (
    fep_map_id          VARCHAR(64) PRIMARY KEY,
    genotype            VARCHAR(64) NOT NULL,
    mue_vs_experiment   DOUBLE,
    n_anchor_pairs      INT NOT NULL DEFAULT 0,
    qc_pass             TINYINT(1),
    rejected_reason     TEXT,                              -- qc_pass=0时必须填
    computed_at         DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT chk_fep_maps_pass CHECK (qc_pass IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE fep_perturbations (
    id                          INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    fep_map_id                  VARCHAR(64) NOT NULL,
    ligand_a                    VARCHAR(64) NOT NULL,
    ligand_b                    VARCHAR(64) NOT NULL,
    predicted_ddg_kcal_mol      DOUBLE NOT NULL,
    uncertainty_kcal_mol        DOUBLE NOT NULL,
    experimental_ddg_kcal_mol   DOUBLE,                     -- 非NULL表示这对是实验锚点
    CONSTRAINT fk_fep_perturbations_map
        FOREIGN KEY (fep_map_id) REFERENCES fep_maps(fep_map_id),
    CONSTRAINT fk_fep_perturbations_ligand_a
        FOREIGN KEY (ligand_a) REFERENCES compounds(compound_id),
    CONSTRAINT fk_fep_perturbations_ligand_b
        FOREIGN KEY (ligand_b) REFERENCES compounds(compound_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================
-- 环节4.3：共价对接(core/covalent.py CovalentDockingResult)
-- ============================================================
CREATE TABLE covalent_docking (
    compound_id                                VARCHAR(64) NOT NULL,
    ensemble_id                                VARCHAR(64) NOT NULL,
    nucleophile_residue                        VARCHAR(32) NOT NULL,
    warhead_carbon_to_nucleophile_distance_a   DOUBLE,
    attack_angle_degree                        DOUBLE,
    attack_geometry_ok                         TINYINT(1),  -- Bürgi-Dunitz判据
    ki_nm                                      DOUBLE,       -- 可逆亲和力，决定选择性
    kinact_per_second                          DOUBLE,       -- 成键速率
    gsh_half_life_hours                        DOUBLE,       -- 弹头本征反应性
    PRIMARY KEY (compound_id, ensemble_id),
    CONSTRAINT fk_covalent_docking_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT fk_covalent_docking_ensemble
        FOREIGN KEY (ensemble_id) REFERENCES structure_ensemble(ensemble_id),
    CONSTRAINT chk_covalent_docking_geom CHECK (attack_geometry_ok IN (0, 1))
    -- kinact_over_ki 不落库，是只读计算属性(kinact/KI)，查询时现算，避免和源字段不同步
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_covalent_docking_ensemble ON covalent_docking (ensemble_id);

-- ============================================================
-- 环节5：选择性引擎(core/selectivity.py SelectivityResult)
-- ============================================================
CREATE TABLE selectivity_results (
    id                      INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    compound_id             VARCHAR(64) NOT NULL,
    genotype                VARCHAR(64) NOT NULL,          -- 目标突变体，比如 del19_C797S
    ddg_target_kcal_mol     DOUBLE,
    ddg_wt_kcal_mol         DOUBLE,
    dddg_kcal_mol           DOUBLE,
    fold_selectivity_pred   DOUBLE,
    uncertainty_kcal_mol    DOUBLE,
    `method`                  VARCHAR(16) NOT NULL,
    structures_used         JSON,                           -- ensemble_id列表
    computed_at             DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_selectivity_results_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_selectivity_results_method CHECK (method IN ('fep', 'mmgbsa', 'ml', 'docking')),
    UNIQUE KEY uq_selectivity_results (compound_id, genotype, method)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================
-- 环节6：ADMET描述符 + 结构警示(core/admet.py)
-- ============================================================
CREATE TABLE admet_descriptors (
    compound_id                 VARCHAR(64) PRIMARY KEY,
    mw                          DOUBLE,
    clogp                       DOUBLE,
    clogd_approx                DOUBLE,                     -- 近似值，见 core/admet.py 顶部声明
    tpsa                        DOUBLE,
    hbd                         INT,
    hba                         INT,
    fsp3                        DOUBLE,
    has_basic_aliphatic_amine   TINYINT(1),
    pka_proxy                   VARCHAR(32),                -- 子结构代理标记，不是数值pKa
    cns_mpo                     DOUBLE,
    cns_mpo_pass_4_0            TINYINT(1),
    herg_risk_proxy             VARCHAR(16),
    computed_at                 DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_admet_descriptors_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_admet_descriptors_amine CHECK (has_basic_aliphatic_amine IN (0, 1)),
    CONSTRAINT chk_admet_descriptors_mpo_pass CHECK (cns_mpo_pass_4_0 IN (0, 1)),
    CONSTRAINT chk_admet_descriptors_herg CHECK (herg_risk_proxy IN ('low', 'medium', 'high'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE structural_alert_hits (
    id              INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    compound_id     VARCHAR(64) NOT NULL,
    rule_name       VARCHAR(64) NOT NULL,                   -- 对应 structural_alerts.yaml 里的key
    severity        VARCHAR(16) NOT NULL,
    `action`          VARCHAR(16) NOT NULL,
    advice          TEXT,
    hit_at          DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_structural_alert_hits_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_structural_alert_hits_action CHECK (action IN ('PASS', 'OPTIMIZE_NEEDED', 'REJECT'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_structural_alert_hits_compound ON structural_alert_hits (compound_id);

-- ============================================================
-- 环节7：合成可及性(core/synthesis.py SynthesisRecord)
-- ============================================================
CREATE TABLE synthesis_records (
    compound_id                 VARCHAR(64) PRIMARY KEY,
    sa_score                    DOUBLE,
    n_steps                     INT,
    longest_linear_sequence     INT,
    key_bb_cas                  VARCHAR(32),
    bb_lead_time_days           INT,
    estimated_cost_per_20mg     DOUBLE,
    chiral_resolution_needed    TINYINT(1),
    scale_up_risk               TEXT,
    fto_status                  VARCHAR(64),                -- in_claim / near_claim / clear / unknown_no_real_patent_data
    CONSTRAINT fk_synthesis_records_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_synthesis_records_chiral CHECK (chiral_resolution_needed IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================
-- 环节8：多目标决策(core/decide.py 的输出留痕，不存算法本身)
-- ============================================================
CREATE TABLE decision_rounds (
    round_id            VARCHAR(32) PRIMARY KEY,           -- 比如 "2026-W40"
    round_date          DATE NOT NULL,
    batch_size          INT NOT NULL,
    quota_json          JSON,                               -- 当轮用的 exploit/explore/hypothesis/control 配额
    notes               TEXT DEFAULT ''
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE decision_batch_members (
    round_id            VARCHAR(32) NOT NULL,
    compound_id         VARCHAR(64) NOT NULL,
    batch_role          VARCHAR(16) NOT NULL,
    is_pareto_optimal   TINYINT(1),
    desirability_score  DOUBLE,
    hypothesis_label    TEXT,                               -- hypothesis_test角色时，药化指定的具体假设描述
    PRIMARY KEY (round_id, compound_id),
    CONSTRAINT fk_decision_batch_members_round
        FOREIGN KEY (round_id) REFERENCES decision_rounds(round_id),
    CONSTRAINT fk_decision_batch_members_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_decision_batch_members_role
        CHECK (batch_role IN ('exploit', 'explore', 'hypothesis_test', 'control')),
    CONSTRAINT chk_decision_batch_members_pareto CHECK (is_pareto_optimal IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_decision_batch_members_compound ON decision_batch_members (compound_id);

-- ============================================================
-- 环节9：DMTA 回传数据(core/feedback.py WetLabResult)
-- ============================================================
CREATE TABLE wet_lab_results (
    id                      INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    -- 身份
    compound_id             VARCHAR(64) NOT NULL,
    batch_id                VARCHAR(64) NOT NULL,
    smiles_as_made          VARCHAR(2000) NOT NULL,         -- 实际做出来的结构，可能和设计的不同，不可省略
    purity_pct              DOUBLE NOT NULL,
    chirality_confirmed     TINYINT(1) NOT NULL,
    -- 活性
    genotype                VARCHAR(64) NOT NULL,
    assay_format            VARCHAR(16) NOT NULL,
    atp_conc_um             DOUBLE,
    readout                 VARCHAR(32) NOT NULL,
    `value`                   DOUBLE NOT NULL,
    `unit`                    VARCHAR(16) NOT NULL,
    censoring               VARCHAR(16) NOT NULL DEFAULT 'point',
    n_replicates            INT NOT NULL,
    log_sd                  DOUBLE,
    assay_date              DATE NOT NULL,
    plate_id                VARCHAR(32),
    -- 选择性
    wt_value                DOUBLE,
    fold_selectivity        DOUBLE,
    kinase_panel_file       VARCHAR(255),
    -- ADME
    sol_ph74                DOUBLE,
    hlm_clint               DOUBLE,
    rlm_clint               DOUBLE,
    ppb_pct                 DOUBLE,
    caco2_papp              DOUBLE,
    mdr1_er                 DOUBLE,
    herg_ic50_um            DOUBLE,
    gsh_adduct_detected     TINYINT(1),
    -- 体内
    `species`                 VARCHAR(32),
    dose_mg_per_kg          DOUBLE,
    `route`                   VARCHAR(32),
    auc                     DOUBLE,
    cmax                    DOUBLE,
    t_half_hours            DOUBLE,
    f_pct                   DOUBLE,
    kp_uu                   DOUBLE,
    -- 元数据
    operator                VARCHAR(64) NOT NULL,
    protocol_version        VARCHAR(32) NOT NULL,
    qc_status                VARCHAR(16) NOT NULL DEFAULT 'pending',
    comments                 TEXT DEFAULT '',
    created_at               DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_wet_lab_results_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT chk_wet_lab_results_purity CHECK (purity_pct BETWEEN 0 AND 100),
    CONSTRAINT chk_wet_lab_results_chirality CHECK (chirality_confirmed IN (0, 1)),
    CONSTRAINT chk_wet_lab_results_assay_format CHECK (assay_format IN ('enzymatic', 'cellular', 'nanobret')),
    CONSTRAINT chk_wet_lab_results_censoring CHECK (censoring IN ('point', 'less_than', 'greater_than')),
    CONSTRAINT chk_wet_lab_results_n_rep CHECK (n_replicates >= 1),
    CONSTRAINT chk_wet_lab_results_gsh CHECK (gsh_adduct_detected IN (0, 1)),
    CONSTRAINT chk_wet_lab_results_qc CHECK (qc_status IN ('pass', 'fail', 'pending'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE INDEX idx_wet_lab_results_compound ON wet_lab_results (compound_id);
CREATE INDEX idx_wet_lab_results_genotype ON wet_lab_results (genotype);

-- 假阳性归因(环节9.2第2条，标签必须结构化进数据库，不能自由文本——季度要统计)
CREATE TABLE false_positive_attributions (
    id                  INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    compound_id         VARCHAR(64) NOT NULL,
    funnel_score_id     INT UNSIGNED,                        -- 具体哪一层的预测被判定为假阳性
    attribution_label   VARCHAR(64) NOT NULL,
    attributed_by        VARCHAR(64) NOT NULL,
    attributed_at        DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    notes                TEXT DEFAULT '',
    CONSTRAINT fk_fp_attr_compound
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
    CONSTRAINT fk_fp_attr_funnel_score
        FOREIGN KEY (funnel_score_id) REFERENCES funnel_scores(id),
    CONSTRAINT chk_fp_attr_label CHECK (attribution_label IN (
        'pose_error', 'protonation_state_error', 'conformer_selection_error',
        'desolvation_underestimated', 'cell_permeability_issue',
        'metabolic_instability', 'compound_degradation_or_purity'
    ))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- 规则治理：precision追踪(环节9.3，core/feedback.py rule_precision_tracking的落库)
-- 列名用 precision_value 而不是 precision——PRECISION 是MySQL保留字
CREATE TABLE rule_precision_history (
    id                          INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    rule_name                   VARCHAR(64) NOT NULL,       -- 对应 structural_alerts.yaml 的key
    review_date                 DATE NOT NULL,
    n_hits                      INT NOT NULL,
    n_true_positive             INT NOT NULL,
    precision_value              DOUBLE NOT NULL,
    n_independent_compounds     INT NOT NULL,
    eligible_for_reject_mode    TINYINT(1) NOT NULL,
    current_mode                VARCHAR(16) NOT NULL,
    CONSTRAINT chk_rule_precision_eligible CHECK (eligible_for_reject_mode IN (0, 1)),
    CONSTRAINT chk_rule_precision_mode CHECK (current_mode IN ('warn', 'reject'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- 前瞻验证时间序列(环节9.2第4条)
CREATE TABLE prospective_validation_rounds (
    round_id            VARCHAR(32) PRIMARY KEY,
    round_date          DATE NOT NULL,
    model_version       VARCHAR(32) NOT NULL,
    n_compounds         INT NOT NULL,
    spearman_rho        DOUBLE,
    mue                 DOUBLE,
    notes               TEXT DEFAULT ''
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
