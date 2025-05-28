# Data Directory

此目录用于存放引物设计流程所需的数据文件 / This directory contains data files required for the primer design pipeline.

## 📁 文件说明 / File Description

### assembly_summary.txt
- **用途 / Purpose**: NCBI基因组数据库摘要文件 / NCBI genome database summary file
- **来源 / Source**: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
- **大小 / Size**: 通常几GB / Usually several GB
- **更新频率 / Update**: 建议定期更新以获取最新基因组信息 / Recommend regular updates for latest genome information

## 📥 如何获取 assembly_summary.txt / How to Get assembly_summary.txt

### 方法1: 直接下载 / Method 1: Direct Download
```bash
# 下载最新的assembly_summary.txt文件
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt -O assembly_summary.txt

# 或使用curl
curl -o assembly_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
```

### 方法2: 使用现有文件 / Method 2: Use Existing File
如果您已经有 `assembly_summary.txt` 文件，请将其复制到此目录：
If you already have an `assembly_summary.txt` file, copy it to this directory:

```bash
cp /path/to/your/assembly_summary.txt src/microbial_primer_design/data/
```

## ⚙️ 配置说明 / Configuration

引物设计流程会自动查找以下位置的 `assembly_summary.txt` 文件：
The primer design pipeline will automatically look for `assembly_summary.txt` in the following locations:

1. **当前工作目录 / Current working directory**
2. **项目data目录 / Project data directory**: `src/microbial_primer_design/data/`
3. **自定义路径 / Custom path**: 通过 `--assembly-summary` 参数指定

## 🔧 使用示例 / Usage Examples

```bash
# 使用data目录中的assembly_summary.txt
cd src/microbial_primer_design
python run_primer_design.py Terrisporobacter --outgroup Intestinibacter

# 使用自定义路径的assembly_summary.txt
python run_primer_design.py Terrisporobacter --assembly-summary /path/to/assembly_summary.txt
```

## ⚠️ 注意事项 / Important Notes

1. **文件大小 / File Size**: `assembly_summary.txt` 文件通常很大（几GB），请确保有足够的磁盘空间
2. **网络下载 / Network Download**: 首次下载可能需要较长时间，建议使用稳定的网络连接
3. **定期更新 / Regular Updates**: 建议定期更新此文件以获取最新的基因组信息
4. **版本控制 / Version Control**: 由于文件很大，建议不要将其加入Git版本控制 