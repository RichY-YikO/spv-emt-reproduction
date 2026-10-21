# SSW event-date provenance

## Date calculation

The SSW event dates were calculated by Shuo Yang using the identification method of Charlton and Polvani (2007; https://doi.org/10.1175/JCLI3996.1). This provenance statement is confirmed by the author; this documentation update does not rerun event detection.

The study manuscript states the basic criterion: reversal of zonal-mean zonal wind at 10 hPa and 60 degrees N. The existing study catalogue contains 47 central dates in 39 July-to-June SSW winters; the other 36 winters are normal winters. The dates in the supplied derived CSVs and the source workbook are retained as calculated by the author.

The supplied reproduction workflow reads the existing event catalogue and validates its dates/workbook checksum; it does not newly derive SSW events from a zonal-wind time series. The author's original event-detection calculation is distinct from that validation step. This update documents author-confirmed provenance; it does not independently verify every detection implementation detail.

## Method reference

Charlton, A. J., and L. M. Polvani (2007), A New Look at Stratospheric Sudden Warmings. Part I: Climatology and Modeling Benchmarks, Journal of Climate, 20(3), 449-469, https://doi.org/10.1175/JCLI3996.1.

An associated methodological correction is available at https://doi.org/10.1175/JCLI-D-11-00348.1. It clarifies the consecutive-westerly-day separation rule. The existence of this reference is not a claim that the original author calculation has been re-audited against every corrected implementation detail.

## 中文说明

具体SSW事件日期由Shuo Yang依据Charlton和Polvani（2007）的方法自行计算。原有47个事件日期、39个SSW冬季与36个normal冬季的分类保持不变。本次仅更新出处说明，不重新识别事件。
