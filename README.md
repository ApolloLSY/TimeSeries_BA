# TimeSeries_BA
start on: 2024.04.12

起因是看到一篇机器学习相关的sci论文把DNA序列当成时间序列数据（一个位点代表一个时间点，作者解释为进化遗传上的时间序列），让我感到非常惊讶可以这样做。虽然这篇论文的设计还是不太能让我认同，不过我忽然开始考虑我在做的课题是否可以当做时间序列数据。

衰老研究领域有一个被拒绝的概念是，人的年龄等于实际经历过的岁月（常见的反例如存活同样时长的人皮肤、器官可以呈现不同的衰老程度，或是患有早衰症的孩子存活10年在生物意义上却可以老得像存活了80年的老人）。研究者们开发了很多种测量衰老的算法，并使用生物年龄这个概念。一大类算法是基于血液检查标志物的组合，另一大类算法是基于DNA甲基化组数据（DNA链上各个部分有不同程度的甲基化，我们通过测量给出每个位点甲基化的概率，并用大量位点数据训练模型）

是否可以把一个人的生物年龄（以DNA甲基化组方法计算得，Horvath1，Hannum2, Lin3, Ying_caus4）当做一种模拟的时间。如果以1岁为间隔，也就是把一个人在生物年龄1,2,3,4,5，……，120，……时（这里不以0起始是有特殊的原因0：在ground zero模型中认为，最初的胚胎和出生时人的生物年龄都不是0，而是遵循一个光滑的近似V形曲线，在接近胚胎中期时出现最低点。由于胚胎中期的生物年龄不易操作，因此本研究设定biological age based pseudo time series的起始点为t=1。），某个位点的甲基化概率当作时间序列来分析。由于一套450k的甲基化组有45万多个观测位点，如果考虑所有位点，它会是一个面板数据。这样建模的意义或许是，可以研究人类DNA甲基化变动是遵循什么样的随机机制，一个人的DNA甲基化组是不是时间依赖的。
还有一个操作层面的问题是，此前的课题是用机器学习处理大量人群横断面DNA甲基化数据，但没有单个人追踪多年的数据。是否可以把同样生物年龄（比如20岁）的人的数据做平均，让n个人代表一个“20岁平均人”，以此类推到每个整岁，构建一个统计意义上的“平均人”从1到120的甲基化面板。如果中心极限定理在这里是适用的，这样应该会近似于一个具体人类在各个岁数的甲基化水平。

GSE41169,GSE40279,GSE51057,GSE42861,GSE51032,GSE73103,GSE69270,GSE64495,GSE30870,GSE52588

0 Gladyshev VN. The Ground Zero of Organismal Life and Aging. Trends Mol Med. 2021 Jan;27(1):11-19. doi: 10.1016/j.molmed.2020.08.012. Epub 2020 Sep 23. PMID: 32980264; PMCID: PMC9183202.(10.1016/j.molmed.2020.08.012)

1 Horvath, S. & Raj, K. DNA methylation-based biomarkers and the epigenetic clock theory of ageing. Nat Rev Genet 19, 371–384 (2018).

2 Hannum, G. et al. Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Molecular Cell 49, 359–367 (2013).

3 Lin, Q. et al. DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy. Aging 8, 394–401 (2016).

4 Ying, K., Liu, H., Tarkhov, A.E. et al. Causality-enriched epigenetic age uncouples damage and adaptation. Nat Aging 4, 231–246 (2024).
