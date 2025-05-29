# 微生物引物设计工具包 / Microbial Primer Design Toolkit

一个用于微生物特异性引物设计的完整流水线工具包。

A comprehensive pipeline toolkit for microbial-specific primer design.

## 🚀 特性 / Features

- **完整的9步自动化流水线** / Complete 9-step automated pipeline
- **并行处理优化** / Parallel processing optimization (89% of workflows)
- **智能质量评分** / Intelligent quality scoring (10-dimensional scoring system)
- **支持NCBI基因组数据** / Support for NCBI genome data
- **一键式使用** / One-command usage

## 📋 系统要求 / System Requirements

- Python 3.7+
- 至少8GB内存 / At least 8GB RAM
- 多核CPU（推荐） / Multi-core CPU (recommended)

## 🔧 安装 / Installation

### 方法1: 使用conda（推荐） / Method 1: Using conda (recommended)

```bash
# 克隆仓库 / Clone repository
git clone https://github.com/YOUR_USERNAME/microbial-primer-design.git
cd microbial-primer-design

# 创建conda环境 / Create conda environment
conda env create -f environment.yml
conda activate primer_design
```

### 方法2: 使用pip / Method 2: Using pip

```bash
# 克隆仓库 / Clone repository
git clone https://github.com/YOUR_USERNAME/microbial-primer-design.git
cd microbial-primer-design

# 安装依赖 / Install dependencies
pip install -r requirements.txt

# 手动安装外部工具 / Manually install external tools
# prokka, roary, mafft, primer3
```

### 外部工具安装 / External Tools Installation

```bash
# 使用conda安装生物信息学工具 / Install bioinformatics tools using conda
conda install -c bioconda prokka roary mafft primer3
```

## 🚀 快速开始 / Quick Start

### 基本使用 / Basic Usage

```bash
# 运行完整流水线 / Run complete pipeline
python run_primer_design.py Terrisporobacter --threads 8

# 包含外群分析 / Include outgroup analysis
python run_primer_design.py Terrisporobacter --outgroup Intestinibacter Clostridium --threads 8

# 限制处理基因数量（用于测试） / Limit gene number (for testing)
python run_primer_design.py Terrisporobacter --max-genes 20 --threads 4
```

### 高级选项 / Advanced Options

```bash
# 从特定步骤开始（断点续传） / Start from specific step (resume)
python run_primer_design.py Terrisporobacter --start-step 3 --threads 8

# 使用高质量模式 / Use high quality mode
python run_primer_design.py Terrisporobacter --high-quality --threads 16

# 自定义assembly_summary.txt路径 / Custom assembly_summary.txt path
python run_primer_design.py Terrisporobacter --assembly-summary /path/to/assembly_summary.txt
```

## 📊 流水线步骤 / Pipeline Steps

1. **📥 基因组下载** / Genome Download - 从NCBI下载基因组 / Download genomes from NCBI
2. **🧬 Prokka注释** / Prokka Annotation - 基因组注释 / Genome annotation
3. **🔄 Roary分析** / Roary Analysis - 泛基因组分析 / Pan-genome analysis
4. **🎯 基因筛选** / Gene Filtering - 筛选特异性基因 / Filter specific genes
5. **📝 序列提取** / Sequence Extraction - 提取基因序列 / Extract gene sequences
6. **🧮 多序列比对** / Multiple Alignment - MAFFT比对 / MAFFT alignment
7. **🧪 引物设计** / Primer Design - Primer3设计 / Primer3 design
8. **📋 结果解析** / Result Parsing - 解析引物信息 / Parse primer information
9. **⭐ 质量评估** / Quality Assessment - 引物质量评分 / Primer quality scoring

## 📁 项目结构 / Project Structure

```
microbial_primer_design/
├── README.md                    # 项目文档 / Project documentation
├── requirements.txt             # Python依赖 / Python dependencies
├── environment.yml              # Conda环境 / Conda environment
├── LICENSE                      # GPL v3许可证 / GPL v3 license
├── run_primer_design.py         # 主运行脚本 / Main script
├── primer_design_toolkit/       # 核心工具包 / Core toolkit
│   ├── primer_pipeline.py       # 主流水线 / Main pipeline
│   ├── genome_downloader.py     # 基因组下载 / Genome downloader
│   ├── primer3_parser.py        # Primer3解析 / Primer3 parser
│   └── quality_ranker.py        # 质量评分 / Quality ranker
└── utils/                       # 独立工具 / Standalone tools
    ├── download_genomes.py      # 基因组下载工具 / Genome download tool
    ├── parse_primer3.py         # Primer3解析工具 / Primer3 parser tool
    └── rank_primers.py          # 引物排序工具 / Primer ranking tool
```

## 📤 输出结果 / Output Results

### 主要输出文件 / Main Output Files

- `final_ranked_primers.csv` - 最终引物结果 / Final primer results
- `genome_summary.csv` - 基因组信息汇总 / Genome information summary
- `specific_genes_info.csv` - 特异性基因信息 / Specific gene information
- `quality_scores.csv` - 质量评分详情 / Quality score details

### 中间文件 / Intermediate Files

- `prokka_results/` - 基因组注释结果 / Genome annotation results
- `roary_results/` - 泛基因组分析结果 / Pan-genome analysis results
- `gene_sequences/` - 基因序列文件 / Gene sequence files
- `alignments/` - 序列比对结果 / Sequence alignment results

## 🔍 独立工具使用 / Standalone Tools Usage

```bash
# 仅下载基因组 / Download genomes only
python utils/download_genomes.py Terrisporobacter --threads 8

# 仅解析Primer3结果 / Parse Primer3 results only
python utils/parse_primer3.py primer3_output.txt parsed_primers.csv

# 仅进行引物质量评分 / Rank primer quality only
python utils/rank_primers.py parsed_primers.csv ranked_primers.csv --summary
```

## 📈 性能优化 / Performance Optimization

- **并行化处理** / Parallel Processing: 89%的流程支持并行 / 89% of workflows support parallelization
- **内存优化** / Memory Optimization: 智能内存管理 / Intelligent memory management
- **速度提升** / Speed Improvement: 相比串行版本提速3-4倍 / 3-4x faster than serial version
- **CPU利用** / CPU Utilization: 充分利用多核CPU / Full utilization of multi-core CPU

## 🔧 故障排除 / Troubleshooting

### 常见问题 / Common Issues

1. **内存不足** / Out of Memory
   ```bash
   # 减少并行线程数 / Reduce parallel threads
   python run_primer_design.py Genus --threads 4
   
   # 限制处理基因数量 / Limit gene number
   python run_primer_design.py Genus --max-genes 50
   ```

2. **外部工具未安装** / External Tools Not Installed
   ```bash
   # 检查工具是否可用 / Check if tools are available
   prokka --version
   roary --version
   mafft --version
   primer3_core --version
   ```

3. **assembly_summary.txt文件问题** / assembly_summary.txt Issues
   ```bash
   # 下载最新的assembly_summary.txt / Download latest assembly_summary.txt
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
   ```

## 📄 许可证 / License

本项目采用GPL v3许可证。详见 [LICENSE](LICENSE) 文件。

This project is licensed under GPL v3. See [LICENSE](LICENSE) file for details.

## 🤝 贡献 / Contributing

欢迎提交问题和拉取请求！

Issues and pull requests are welcome!

## 📞 联系 / Contact

如有问题，请提交GitHub Issue。

For questions, please submit a GitHub Issue.
