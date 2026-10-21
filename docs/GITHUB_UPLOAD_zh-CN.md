# GitHub登录与上传步骤

## 是否需要登录

整理和检查本地发布包不需要账号。上传到自己的GitHub仓库需要GitHub账号，并在GitHub Desktop或命令行中完成认证。登录在GitHub官方页面或客户端完成，不需要把密码、令牌或验证码发给别人。

本次只交付本地文件夹和ZIP，没有创建远程仓库、提交代码或上传任何内容。

## 推荐：GitHub Desktop

1. 完成`RELEASE_CHECKLIST.md`中的署名、许可证和数据来源确认。若仅先做协作检查，可先建立私有仓库，公开时再核对清单。
2. 安装GitHub Desktop并使用你的GitHub账号登录。
3. 解压交付ZIP，选择其中的`spv-emt-reproduction-v1.0.0`文件夹。不要将ZIP本身当作仓库的唯一文件上传。
4. 在GitHub Desktop选择 **File → Add local repository** 并选择该文件夹。若提示尚不是Git仓库，点击提示中的 **create a repository here**，在该目录创建本地仓库；不要覆盖现有README、gitignore或正式LICENSE，不另选自动许可证。
5. 查看Changes列表，确认包含代码、缓存、文档、最终图，排除了`outputs/`、原始`data/raw/`、Python缓存、备份、账号文件和交付ZIP。`.gitignore`已提供这些规则。
6. 若还未创建初始提交，输入`Prepare reproducible SPV EMT release`并点击 **Commit to main**。
7. 点击 **Publish repository**，仓库名用`spv-emt-reproduction`或你确认的名称。是否保留 **Keep this code private** 由实际审稿要求决定；取消勾选即公开。点击 **Publish Repository**完成上传。
8. 在浏览器打开仓库，确认README正常显示、文件完整、许可和引用信息正确。后续修改使用Commit后再Push origin。

GitHub Desktop的发布步骤参考：[官方说明](https://docs.github.com/en/desktop/adding-and-cloning-repositories/adding-an-existing-project-to-github-using-github-desktop)。

## 为什么不建议网页拖拽整个包

GitHub网页直接上传每个文件限制25 MiB；Git方式可提交至100 MiB的单文件，超过100 MiB需要Git LFS。本包的最终Fig.1缓存约34.7 MiB，超过网页限额但低于Git单文件上限，因此使用GitHub Desktop或Git即可，无需仅为当前文件启用LFS。整个ZIP大小不等于Git单文件大小；不要把大ZIP提交进仓库。

限制和上传方式参考：[GitHub官方文件上传说明](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)。

## 可选：Git命令行

在发布包根目录打开终端，先通过官方方式完成GitHub认证，再执行：

```bash
git init -b main
git add .
git status
git commit -m "Prepare reproducible SPV EMT release"
```

如果Git要求姓名和邮箱，请填写你希望显示在提交记录中的署名；它不等于论文软件作者名单。然后在GitHub网页创建一个空仓库，不自动添加README或LICENSE。将下面占位地址换成你的真实仓库地址：

```bash
git remote add origin https://github.com/YOUR_ACCOUNT/spv-emt-reproduction.git
git push -u origin main
```

仓库发布后可创建一个固定版本Release。若审稿要求可长期引用的软件存档，可再按所选归档平台流程生成DOI，并回填正式CITATION.cff及论文数据/代码可用性声明；当前没有生成任何DOI。
