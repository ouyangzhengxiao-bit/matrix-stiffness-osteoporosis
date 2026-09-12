# eBMD 遗传学外部验证：预设分析方案 v1

**冻结日期：2026-09-08。** 本方案在将刚度元分析基因与 eBMD 基因统计量合并之前冻结。

## 数据与主要问题

- 遗传数据：2026 年 Nature Genetics 研究的 UK Biobank eBMD MAGMA 全基因组基因检验（448,010 人，Supplementary Table 7e）。
- 暴露侧数据：此前 7 个独立刚度转录组研究的随机效应元分析。
- 主要问题：元分析预先固定的 200 基因刚度 signature（正向 100、负向 100 的并集）是否比可比背景基因具有更强的 eBMD 基因关联。
- MAGMA `ZSTAT` 只表示关联强度，不提供可与表达方向对应的遗传效应方向，因此主要分析不声称“升高/降低骨密度”。

## 主要检验

主要统计量为 signature 内基因的平均 MAGMA `ZSTAT`。背景抽样同时匹配 gene length 十分位和 `NSNPS` 十分位；从非 signature 蛋白编码基因中无放回抽取 20,000 次。方向性备择假设为 signature 平均 ZSTAT 更高，报告单侧经验 P 值和标准化效应。

主要检验只有一个，不与探索性分析合并校正。不能用显著 eBMD 基因阈值反向筛选或修改刚度 signature。

## 预设次要分析

1. 正向与负向 signature 分别做匹配竞争性检验，两个 P 值 BH-FDR。
2. 刚度元分析 `abs(pooled_effect)` 与 MAGMA ZSTAT 的 Spearman 相关。
3. 作者定义的 5 个人 BMSC 单细胞富集 marker 集合分别做匹配竞争性检验，5 个 P 值 BH-FDR。
4. 刚度 signature 与显著 eBMD 基因（论文阈值 `P < 2.5e-6`）的超几何重叠，仅作为可解释性分析。
5. 回归敏感性分析：`ZSTAT ~ signature + log(gene length) + log(NSNPS) + log(NPARAM) + chromosome`；报告 signature 系数与 HC3 稳健标准误。
6. ATAC promoter delta 与 MAGMA ZSTAT 的相关及三层联合优先级仅作探索性分析。

## 判定与表述

- 主要检验 `P < 0.05` 且标准化效应为正：支持刚度 signature 富集 eBMD 遗传关联。
- 方向为正但未达显著：趋势证据。
- 方向不为正：不支持。
- 不论结果如何，观察性整合均不表述为因果机制或临床生物标志物。
