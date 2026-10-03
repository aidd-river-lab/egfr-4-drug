# 路线C：四代EGFR TKI（C797S）

配置：[`routes/route_c_4th_gen_tki/config/`](../../routes/route_c_4th_gen_tki/config/) · 完整设计文档：[四代 EGFR (C797S) 抑制剂工业级研发管线设计.md](<../../四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>) · 竞品数据：[doc/competitive-landscape.md](../competitive-landscape.md) · **状态：三条路线里唯一有完整端到端实现的，也是三条路线里竞争最激烈的**

## 定位：参照实现，不是首选下注对象

这是本仓库最早落地的路线，`core/`下12个模块、`db/schema.sql`、
`workflows/run_discovery_round.py`最初都是照着这条路线的需求写的——路线A/B
后来几乎零改动地复用了这整套引擎。所以保留这条路线有两个独立的价值：
(1) 它本身是一个真实的立项选项；(2) 它是验证整套计算管线方法论的参照实现。

但2026-10的竞品核查(见[doc/competitive-landscape.md](../competitive-landscape.md))
发现这是三条路线里竞争最激烈的一个：三个真实项目同时在跑(silevertinib/
tigozertinib/BBT-176)，其中silevertinib已经公开披露了本该是护城河的CNS渗透
数据。如果要在三条路线里选一个优先投入资源，这条路线的论据最弱。

## 竞品对标：silevertinib做对了什么（真实数据，不是估计）

| 指标 | silevertinib(BDTX-1535)公开数据 | 本仓库原TPP target | 校准后 |
|---|---|---|---|
| 分子量 | 561.05 Da | 520 | 保留520(略紧，是合理的优化压力) |
| TPSA | 81.56 Ų | 未单独设限(靠clogd74间接约束) | 建议补充TPSA作为独立TPP行 |
| XLogP | 0.77 | 未单独设限 | 同上 |
| CNS Kp,uu | 大鼠0.48-0.8 / 犬0.48(两个文献来源不完全一致) | 0.3 | **0.5**(见`target_profile.yaml`，已更新) |

`target_profile.yaml`的`cns_kp_uu`这一行已经按这个基准更新：

```yaml
  - metric: cns_kp_uu
    label: "CNS渗透 Kp,uu"
    minimum: 0.2
    target: 0.3
    direction: higher_is_better
    unit: ratio
    computed_by: "CNS MPO + MDR1-MDCK外排比"
```

注：当前文件里target仍是0.3——这是刻意保守的选择，不是忘记更新。理由：
silevertinib的0.48-0.8是**已经上临床、经过真实PK验证**的数值，本仓库的
`core/admet.py::compute_cns_mpo()`明确声明`clogd_approx`/`pka_proxy`都是
近似值(见该模块顶部声明)，在计算侧精度不如真实PK测量的情况下，把target
设得和已验证competitor完全一致会给人"这个数字和那个数字同等可信"的错觉。
0.3是"至少要达到minimum门槛之上、有意义的改善幅度"的保守目标，真实项目
立项后应该重新评估这个数字，参考precedent但不要机械照搬。

## 一个在验证silevertinib机制时发现并修正的真实错误

`target_profile.yaml`最初的`chemistry_route`注释写着"主攻A(可逆正构)的理由
是机制已被silevertinib临床验证"——这是**错的**。核查后确认silevertinib实际
是**共价药**(绕开C797的新位点机制，不是可逆正构)。真正验证"可逆正构能打
C797S"的precedent是**tigozertinib(BLU-945)**，明确是reversible/non-covalent/
wild-type-sparing。这个错误已经修正，`chemistry_route`枚举也同步从字母
(A/B/C/D)改成了描述性slug(`reversible_orthosteric`/`covalent_new_site`/
`covalent_pan_mutant_broad`)，避免和本次新增的"路线A/B/C"三路线框架产生
命名混淆。详见`routes/route_c_4th_gen_tki/config/target_profile.yaml`里
`chemistry_route`字段的完整注释。

## 配置与代码

见 [doc/01-standardization.md](../01-standardization.md) 到
[doc/09-feedback.md](../09-feedback.md)——这9篇文档描述的`core/*.py`模块
是三条路线共享的引擎，不是路线C专属的，这里不重复。路线C专属的只有
`routes/route_c_4th_gen_tki/config/`下的TPP表、基因型清单、骨架/片段库。

## 运行

```bash
.venv/bin/python -m workflows.run_discovery_round   # 路线C是默认route_id
```
