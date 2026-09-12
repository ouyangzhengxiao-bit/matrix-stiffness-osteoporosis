# 人骨单细胞与空间组织定位：预设方案 v1

**冻结日期：2026-09-08。** 在合并刚度元分析结果与 Nature Genetics 2026 人骨基因程序前冻结。

## 数据与假设

- 人股骨 scRNA-seq：125,063 个细胞、20 个作者定义细胞群的差异表达基因程序（Supplementary Table 2c）。
- 独立空间转录组：11 个作者定义细胞群的基因程序（Supplementary Table 2d）。
- 主要假设：刚度反应强度 `abs(meta z)` 集中在人骨 scRNA 的成骨谱系程序，包括两个 MSC 群、pre-osteoblasts 和 mature osteoblasts。

## 主要检验

将四个成骨谱系群的作者定义程序取并集。主要统计量为该并集基因的平均 `abs(meta z)`。从非该并集背景中，按元分析纳入研究数 `k` 和 `se_hk` 十分位匹配抽样 20,000 次，做单侧竞争性检验。

## 次要与复制分析

1. 对 20 个 scRNA 细胞群分别检验 signed meta z 与 abs(meta z)，每类 20 个 P 值分别 BH-FDR。
2. 在独立空间数据中，将 MSC、Osteo-MSC 和 Osteoblasts 程序取并集，以相同方法检验 abs(meta z)，作为预设复制。
3. 对 11 个空间细胞群分别检验 signed/absolute meta z，并分别 BH-FDR。
4. 用程序的连续 average log2FC 与 meta z/abs(meta z) 做群内 Spearman 相关，作为探索性结果。

结果只用于组织/细胞状态定位，不把 marker 富集解释为细胞比例变化或因果机制。
