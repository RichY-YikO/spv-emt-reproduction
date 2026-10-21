# SPV本征微观态复现发布包

这是2026-10-19准备的本地归档候选包，版本为`v1.0.0`，不是已经发布的记录。软件作者和版权持有人均为Shuo Yang；没有登录、创建或上传任何Zenodo/GitHub记录。

归档标题：Code and reproduction materials for “State-resolved multimodal contributions to stratospheric polar vortex predictability”。正式发布日期待实际发布时填写。

## 文件位置

- `code/`：最终计算、绘图和验证代码，以及必要的NPZ缓存和CSV。
- `figures/pdf/`：14幅最终矢量PDF，逐字节保留当前最终版。
- `figures/png/`：14幅600 dpi预览图。
- `docs/`：方法、数据要求、图件对应关系、上传步骤和发布前待补信息。
- `validation/`：文件校验清单、代码清理记录和本次实测核验报告。

解释性#注释与模块、函数、类的docstring均已删除，颜色字符串中的#不是注释。方法说明移至文档，命令行帮助仍可使用。旧自定义相关统计版本和历史备份不在发布包中；必要的原始S模态参考缓存保留用于验证，不作为最终Fig.1数据来源。

## 安装与重绘

在本发布包根目录打开终端，建议使用Python 3.10：

```powershell
python -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
python verify_release.py
python reproduce.py
```

也可用`conda env create -f environment.yml`创建环境，再执行`conda activate spv-emt`。

```powershell
python reproduce.py --figures Fig1 FigS1
python reproduce.py --figures Fig3 Fig4 Fig5 --output-dir my-output
python code/redraw_final.py --figures Fig2
python reproduce.py --help
```

默认输出到`outputs/pdf/`和`outputs/png/`。不会覆盖随包提供的最终图，也不会重算本征分解、预测或999次替代检验。未指定图号时重绘14幅图；`Fig1a`可单独绘制Fig.1上方面板。

## 科学定义

最终Fig.1、Fig.S1采用标准Pearson相关：取反原始时间系数后，减去取反后时间序列在该滞后有效配对样本内的自身均值；SPV也在同一配对样本中心化。绘图时不额外翻转。六幅空间图使用原来的`-U_raw`，不改变方向。

Fig.1的灰带考虑0–60天内最大绝对相关的滞后选择；Fig.S1的灰带是固定滞后的零假设包络，二者不混用。均使用999次替代与既定随机种子。

## 验证

```powershell
python verify_release.py --mode numerical
python verify_release.py --mode si
python -m pip install -r requirements-validation.txt
python verify_release.py --mode pdf
```

数值核验使用缓存及独立相关公式，不读取原始NetCDF。`--mode all`依次运行全部验证。完整原始资料重算步骤与文件格式见`docs/DATA_AND_REPRODUCTION.md`；本次不执行耗时的全流程重算。

地图第一次绘制可能需要下载Cartopy海岸线数据。Times New Roman字体需自行安装，发布包不附带字体。不同环境可能改变字体渲染或重新分解时模态的任意正负号；随包缓存固定了已确认的最终方向。

## 来源说明与后续归档

`LICENSE`已填写Copyright (c) 2026 Shuo Yang，适用于原始软件代码；`CITATION.cff`已列软件作者Shuo Yang和版本v1.0.0。ERA5来源已根据作者确认补齐，DOI为https://doi.org/10.24381/cds.bd0915c6。SSW事件日期由Shuo Yang依据Charlton与Polvani（2007）的方法自行计算；数据、图件与第三方材料的再分发许可不由软件MIT许可自动决定。

整理本地文件不需要GitHub账号；上传到自己的仓库需要账号登录或认证。按`docs/GITHUB_UPLOAD_zh-CN.md`操作即可，无需向他人提供密码或验证码。GitHub发布后，如需论文引用的固定版本，可再创建Release，并另行考虑归档DOI；本包没有编造DOI。Zenodo可直接上传ZIP，无需先上传GitHub，填写说明见`docs/ZENODO_UPLOAD_zh-CN.md`；本次只准备本地压缩包。

2026-10-19排版更新：Fig.1六幅空间图尺寸已与Fig.2统一（尺寸参数1.56英寸），仍为两排三列；相关值、空间模态数值及其他图件保持不变。

## 2026-10-20目录精简

`code/fig1_fig2_results/`只保留最终S/T缓存和三份CSV；最终S缓存为`S_final_pearson.npz`，`S_top100.csv`对应最终Fig.1。原定向S缓存及SPV日序列移到`code/reference_inputs/`，供Fig.4重新计算和数值验证使用。两个无读取用途的气候态缓存及旧S相关CSV已从发布目录删除。

Table S1 Excel在`figures/TableS1.xlsx`；全部24个数值已与CSV核对。图件排版JSON在`validation/layout/`，Fig.1计算说明在`validation/fig1_pearson/`，并非科学计算缓存。重绘时对应JSON写到`outputs/validation/layout/`。历史验证记录保留其原执行信息，当前路径以新说明和校验清单为准。
