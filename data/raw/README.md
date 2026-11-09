# 原始 GBD CSV：只放在本地

把以下两份文件存入本目录，不要更改文件内容：
- `IHME-GBD_2023_DATA-b40ffc46-1.csv`：GBD 2023，Level 2，19 类；1980–2023 年，六个五岁年龄组，男女，Deaths Number 和 Rate。
- `IHME-GBD_2023_DATA-97ee48e8-1.csv`：GBD 2023，Level 1，3 大类；其余条件一致。

两份文件都来自 IHME [GBD Results Tool](https://vizhub.healthdata.org/gbd-results/)。

核验结果：全部 528 个年龄 × 年份 × 性别单元上，Level 1 和 Level 2 死亡人数合计一致；不需要另外下载 `All causes` 文件。

原始 CSV 在 `.gitignore` 中排除，不公开上传。R 脚本可以再生成中间结果及原文相应图表。
