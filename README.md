<<<<<<< HEAD
# 微生物引物设计工具包 / Microbial Primer Design Toolkit

一个用于微生物特异性引物设计的完整流水线工具包。

A comprehensive pipeline toolkit for microbial-specific primer design.

## 🚀 快速开始 / Quick Start

### 系统要求 / System Requirements

- Python 3.7+
- conda 或 pip
- 至少8GB内存
- 多核CPU（推荐）

### 安装 / Installation

#### 方法1: 使用conda（推荐）
```bash
# 克隆仓库
git clone https://github.com/YOUR_USERNAME/microbial-primer-design.git
cd microbial-primer-design

# 创建conda环境
conda env create -f environment.yml
conda activate primer_design
```

#### 方法2: 使用pip
```bash
# 克隆仓库
git clone https://github.com/YOUR_USERNAME/microbial-primer-design.git
cd microbial-primer-design

# 安装依赖
pip install -r requirements.txt

# 手动安装外部工具
# prokka, roary, mafft, primer3
```

### 基本使用 / Basic Usage

```bash
# 运行完整流水线
python run_primer_design.py --genus Terrisporobacter --threads 8

# 从特定步骤开始
python run_primer_design.py --genus Terrisporobacter --start-step 3

# 设置最大基因数量
python run_primer_design.py --genus Terrisporobacter --max-genes 20

# 查看帮助
python run_primer_design.py --help
```

## 📁 项目结构 / Project Structure

```
microbial_primer_design/
├── README.md                    # 项目文档
├── requirements.txt             # Python依赖
├── environment.yml              # Conda环境
├── LICENSE                      # MIT许可证
├── .gitignore                   # Git忽略规则
├── run_primer_design.py         # 主运行脚本
├── primer_design_toolkit/       # 核心工具包
│   ├── __init__.py
│   ├── primer_pipeline.py       # 主流水线
│   ├── genome_downloader.py     # 基因组下载
│   ├── primer3_parser.py        # Primer3解析
│   └── quality_ranker.py        # 质量评分
├── utils/                       # 独立工具
│   ├── parse_primer3.py         # 解析Primer3结果
│   ├── rank_primers.py          # 引物质量排序
│   └── download_genomes.py      # 基因组下载工具
└── data/                        # 数据目录
    └── assembly_summary.txt     # NCBI基因组信息（用户提供）
```

## 📊 流水线步骤 / Pipeline Steps

### 完整9步流水线：

1. **📥 基因组下载** - 下载目标和外群基因组
   - 基于NCBI assembly_summary.txt
   - 支持并行下载
   
2. **🧬 Prokka注释** - 基因组注释
   - 并行处理多个基因组
   - 生成GFF和蛋白质序列文件

3. **🔄 Roary分析** - 泛基因组分析  
   - 识别核心基因和特异性基因
   - 高质量参数优化

4. **🎯 基因筛选** - 筛选特异性基因
   - 只在目标菌属中存在
   - 外群中不存在
   - 优先选择hypothetical protein

5. **📝 序列提取** - 提取基因序列
   - 从注释文件提取序列
   - 并行处理

6. **🧮 多序列比对** - MAFFT比对和保守区域识别
   - 高质量序列比对
   - 识别80-400bp保守区域

7. **🧪 引物设计** - Primer3设计
   - 针对保守区域设计引物
   - 并行处理多个区域

8. **📋 结果解析** - 解析引物信息
   - 提取引物序列和参数
   - 生成结构化数据

9. **⭐ 质量评估** - 引物质量评分
   - 10维质量评分体系
   - 智能排序和筛选

## 🔧 高级选项 / Advanced Options

### 命令行参数

```bash
python run_primer_design.py [选项]

必需参数:
  --genus GENUS               目标菌属名称

可选参数:
  --threads THREADS           并行线程数 (默认: 4)
  --start-step STEP          从指定步骤开始 (1-9)
  --max-genes N              最大处理基因数 (默认: 50)
  --working-dir DIR          工作目录
  --assembly-summary PATH    assembly_summary.txt路径
```

### 数据准备

1. **获取assembly_summary.txt**：
   ```bash
   # 从NCBI下载
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
   
   # 放置到data目录
   cp assembly_summary.txt data/
   ```

2. **目录结构**：
   - 确保data目录存在
   - assembly_summary.txt正确放置

## 📈 性能优化 / Performance

- **并行化**: 89%的流程支持并行处理
- **内存优化**: 智能内存管理，支持大规模数据
- **速度提升**: 相比串行版本提速3-4倍
- **CPU利用**: 充分利用多核CPU资源

## 🎯 典型用例 / Use Cases

### 用例1: 新菌属引物设计
```bash
python run_primer_design.py --genus NewGenus --threads 16 --max-genes 100
```

### 用例2: 快速验证
```bash
python run_primer_design.py --genus TestGenus --threads 4 --max-genes 10
```

### 用例3: 从特定步骤恢复
```bash
python run_primer_design.py --genus Terrisporobacter --start-step 6 --threads 8
```

## 📤 输出结果 / Output Results

### 主要输出文件：

- `final_ranked_primers.csv` - 最终引物结果
- `genome_summary.csv` - 基因组信息汇总
- `specific_genes_info.csv` - 特异性基因信息
- `primer3_results/` - Primer3原始结果
- `quality_scores.csv` - 质量评分详情

### 中间文件：

- `prokka_results/` - 基因组注释结果
- `roary_results/` - 泛基因组分析结果
- `gene_sequences/` - 基因序列文件
- `alignments/` - 序列比对结果

## 🔍 问题排查 / Troubleshooting

### 常见问题：

1. **内存不足**: 减少threads数量或max-genes数量
2. **找不到基因组**: 检查assembly_summary.txt是否正确
3. **Prokka失败**: 确保conda环境正确激活
4. **Primer3错误**: 检查保守区域是否合适

## 📄 许可证 / License

MIT License - 详见LICENSE文件

## 🤝 贡献 / Contributing

欢迎提交Issue和Pull Request！

## 📞 联系 / Contact

如有问题，请创建GitHub Issue。

---

**🎯 一键运行，轻松获得高质量微生物特异性引物！**
=======
# 微生物引物设计完整流程 (全面并行化优化版)

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Performance](https://img.shields.io/badge/Performance-3.4x_faster-red.svg)](#性能提升)

一个完整的微生物引物设计流程，从基因组下载到引物质量评分的全自动化流程。通过全面并行化优化

## 🚀 主要特性

### 📋 完整流程
1. **基因组下载** (并行优化) - 自动下载目标和外群基因组
2. **Prokka注释** (并行优化) - 高度并行基因组注释
3. **Roary分析** (部分并行) - 泛基因组分析
4. **基因筛选** (并行优化) - 筛选目标特异性基因
5. **序列提取** (并行优化) - 双层并行序列提取
6. **多序列比对** (并行优化) - 完全并行化比对和保守区域识别
7. **引物设计** (并行优化) - 并行Primer3设计
8. **结果解析** (并行优化) - 并行解析和结构化提取
9. **质量评分** (并行优化) - 多维度并行评分

### 🎯 智能特性
- **断点续传**: 支持从任意步骤开始运行
- **智能基因选择**: 优先选择hypothetical protein基因
- **自动外群识别**: 基于前缀自动分类目标和外群基因组
- **灵活配置**: 支持线程数、基因数量、质量模式等多种配置

## 📦 快速安装

### 方法1: 一键安装 (推荐)
```bash
# 下载项目
git clone <repository_url>
cd src/microbial_primer_design

# 运行自动安装脚本
bash install.sh
```

### 方法2: 手动安装
```bash
# 1. 安装conda环境
conda env create -f environment.yml
conda activate primer_design

# 2. 创建专用环境
conda create -n prokka -c conda-forge -c bioconda prokka
conda create -n roary -c conda-forge -c bioconda roary

# 3. 下载数据文件
mkdir -p data
wget -O data/assembly_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

# 4. 验证安装
python test_installation.py
```

### 方法3: 使用pip安装Python依赖
```bash
pip install -r requirements.txt
```

## 🔧 系统要求

### 硬件要求
- **CPU**: 16核心以上推荐 (32核心最佳)
- **内存**: 32GB以上推荐 (64GB最佳)
- **存储**: 100GB以上可用空间，SSD推荐
- **网络**: 稳定的互联网连接

### 软件依赖
- **Python**: 3.7+ (推荐3.8或3.9)
- **Conda**: Miniconda或Anaconda
- **外部工具**: Prokka, Roary, MAFFT, Primer3

详细安装指南请查看 [INSTALLATION_GUIDE.md](INSTALLATION_GUIDE.md)

## 🚀 快速开始

### 基本用法
```bash
# 激活环境
conda activate primer_design

# 基本运行
python run_primer_design.py Terrisporobacter --outgroup Intestinibacter

# 查看帮助
python run_primer_design.py --help
```

### 推荐配置

#### 大规模数据处理
```bash
python run_primer_design.py Terrisporobacter \
    --outgroup Intestinibacter \
    --threads 50 \
    --max-genes 100 \
    --high-quality
```

#### 快速测试
```bash
python run_primer_design.py Terrisporobacter \
    --outgroup Intestinibacter \
    --threads 20 \
    --max-genes 20 \
    --fast
```

#### 断点续传
```bash
# 从步骤4开始 (跳过下载和注释)
python run_primer_design.py Terrisporobacter \
    --start-step 4 \
    --threads 30
```



## 📁 输出文件

### 主要结果文件
```
Genus_Name/
├── pipeline.log                           # 运行日志
├── Genus_Name_genome_summary.csv          # 基因组摘要信息
├── data/                                  # 下载的基因组文件
├── prokka_results/                        # Prokka注释结果
├── roary_results/                         # Roary泛基因组分析结果
├── specific_genes.txt                     # 特异性基因列表
├── specific_genes_detailed.csv            # 详细基因信息
├── gene_sequences/                        # 提取的基因序列
├── alignments/                            # 多序列比对结果
├── conserved_regions.txt                  # 保守区域信息
└── primer3_results/                       # 引物设计结果
    ├── primer3_output.txt                 # Primer3原始输出
    ├── parsed_primers.csv                 # 解析后的引物信息
    └── ranked_primers.csv                 # 最终排序的引物结果 ⭐
```

### 关键结果文件说明
- **`ranked_primers.csv`**: 最终的引物设计结果，包含质量评分和排序
- **`specific_genes_detailed.csv`**: 筛选出的目标特异性基因详细信息
- **`conserved_regions.txt`**: 识别出的保守区域信息
- **`pipeline.log`**: 完整的运行日志，包含所有步骤的详细信息

## 🎯 参数说明

### 必需参数
- `genus`: 目标genus名称

### 可选参数
- `--outgroup`: 外群genus名称列表
- `--threads`: 并行线程数 (默认: 4, 推荐: 20-50)
- `--max-genes`: 最大处理基因数量 (默认: 全部, 推荐: 20-100)
- `--start-step`: 开始步骤 (1-9, 支持断点续传)
- `--fast`: 快速模式 (默认)
- `--high-quality`: 高质量模式 (更慢但结果更好)
- `--assembly-summary`: assembly_summary.txt文件路径

### 模式选择
- **快速模式** (`--fast`): 使用Roary的`-n`参数，适合初步分析
- **高质量模式** (`--high-quality`): 使用Roary的`--mafft`参数，适合最终分析

## 🧬 基因选择策略

流程会智能选择特异性基因用于引物设计：

1. **优先选择**: hypothetical protein基因 (≥5个时)
2. **备选方案**: 如果hypothetical protein不足，使用所有特异性基因
3. **数量控制**: 通过`--max-genes`参数控制处理的基因数量
4. **外群识别**: 自动识别目标genus前缀开头的为目标基因组，其他为外群

例如：`Terrisporobacter` → `Ter_`开头为目标，其他(如`Int_`)为外群

## 🐛 故障排除

### 常见问题

#### 1. 导入错误
```bash
ModuleNotFoundError: No module named 'primer_design_toolkit.pipeline'
```
**解决方案**: 检查`__init__.py`文件中的导入路径是否正确

#### 2. Primer3失败
```bash
ERROR: ❌ Primer3运行失败: Command '['primer3_core']' returned non-zero exit status 255
```
**解决方案**: 
- 检查primer3是否正确安装: `primer3_core --help`
- 代码已移除有问题的配置路径，使用默认配置

#### 3. 内存不足
**解决方案**:
- 减少线程数: `--threads 10`
- 限制基因数量: `--max-genes 20`
- 监控内存使用: `htop` 或 `free -h`

#### 4. 网络下载失败
**解决方案**:
- 检查网络连接: `ping ftp.ncbi.nlm.nih.gov`
- 手动下载assembly_summary.txt文件
- 使用代理设置 (如果需要)

### 获取帮助
1. 运行 `python test_installation.py` 检查环境
2. 查看 `pipeline.log` 获取详细错误信息
3. 查看 [INSTALLATION_GUIDE.md](INSTALLATION_GUIDE.md) 获取详细安装指南
4. 查看 [PARALLELIZATION_ANALYSIS.md](PARALLELIZATION_ANALYSIS.md) 了解性能优化详情

## 📚 文档

- **[INSTALLATION_GUIDE.md](INSTALLATION_GUIDE.md)**: 详细的安装和依赖指南
- **[PARALLELIZATION_ANALYSIS.md](PARALLELIZATION_ANALYSIS.md)**: 并行化优化分析报告
- **[requirements.txt](requirements.txt)**: Python依赖包列表
- **[environment.yml](environment.yml)**: Conda环境配置文件

## 🔬 示例分析

### 完整示例: Terrisporobacter引物设计
```bash
# 1. 激活环境
conda activate primer_design

# 2. 运行完整流程
python run_primer_design.py Terrisporobacter \
    --outgroup Intestinibacter \
    --threads 30 \
    --max-genes 50 \
    --high-quality

# 3. 查看结果
ls Terrisporobacter/primer3_results/ranked_primers.csv
```

### 预期输出
```
🚀 微生物引物设计流程 (全面并行化优化版)
============================================================
🎯 目标genus: Terrisporobacter
🔗 外群: Intestinibacter
🧵 并行线程数: 30
📊 最大基因数: 50
⚡ 分析模式: 高质量模式 (--mafft)
🔄 开始步骤: 1
============================================================

[执行过程...]

🎉 流程执行成功！
📊 性能提升: 相比串行版本提升3.4-3.9倍
⏱️  总耗时: 约30-80分钟 (取决于数据规模和硬件配置)
```

## 📈 版本信息

- **当前版本**: v2.0 (全面并行化版)
- **Python要求**: 3.7+
- **主要依赖**: pandas≥1.3.0, biopython≥1.79, numpy≥1.21.0
- **外部工具**: prokka 1.14.6, roary 3.13.0, mafft 7.505, primer3 2.6.1

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！

## 📄 许可证

本项目采用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 🙏 致谢

感谢以下开源项目的支持：
- [Prokka](https://github.com/tseemann/prokka) - 基因组注释
- [Roary](https://github.com/sanger-pathogens/Roary) - 泛基因组分析
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) - 多序列比对
- [Primer3](https://github.com/primer3-org/primer3) - 引物设计
- [Biopython](https://biopython.org/) - 生物信息学工具

---

**🚀 享受3.4-3.9倍的性能提升，让引物设计更高效！** 

# 微生物引物设计工具包 / Microbial Primer Design Toolkit

一个用于微生物特异性引物设计的完整流水线工具包，支持从基因组下载到引物质量评估的全流程自动化处理。

A comprehensive pipeline toolkit for microbial-specific primer design, supporting automated processing from genome download to primer quality assessment.

## 🌐 项目语言规范 / Project Language Standards

### 对话与交流 / Communication
- **中文**: 所有技术讨论、问题解答、文档说明均使用中文
- **Chinese**: All technical discussions, Q&A, and documentation explanations in Chinese

### 代码规范 / Code Standards
- **注释 Comments**: 中英双语对照 / Bilingual (Chinese/English)
  ```python
  # 初始化引物评分器 / Initialize primer ranker
  ranker = PrimerQualityRanker(threads=4)
  ```

- **输出信息 Output**: 英文 / English only
  ```python
  print("Primer quality evaluation completed!")
  self.logger.info("Loading primer data from file")
  ```

- **变量命名 Variables**: 英文 / English only
  ```python
  quality_score = 85.5  # ✅ 正确 / Correct
  质量得分 = 85.5       # ❌ 错误 / Incorrect
  ```

### 文档规范 / Documentation Standards
- **README**: 中文为主，关键术语英文标注 / Chinese primary, key terms in English
- **技术文档**: 中英双语对照 / Bilingual technical documentation
- **代码文档**: 中英双语文档字符串 / Bilingual docstrings

## 📁 项目结构 / Project Structure

```
src/microbial_primer_design/
├── .project_config              # 项目配置文件 / Project configuration
├── CODING_STANDARDS.md          # 编码标准文档 / Coding standards documentation
├── README.md                    # 项目说明 / Project documentation
├── run_primer_design.py         # 主运行脚本 / Main execution script
├── primer_design_toolkit/       # 核心工具包 / Core toolkit
│   ├── __init__.py
│   ├── config_manager.py        # 配置管理器 / Configuration manager
│   ├── pipeline.py              # 流水线引擎 / Pipeline engine
│   ├── genome_downloader.py     # 基因组下载器 / Genome downloader
│   ├── primer3_parser.py        # Primer3解析器 / Primer3 parser
│   └── quality_ranker.py        # 质量评分器 / Quality ranker
├── utils/                       # 独立工具集 / Independent utilities
│   ├── parse_primer3.py         # Primer3结果解析 / Primer3 result parsing
│   ├── rank_primers.py          # 引物质量排序 / Primer quality ranking
│   └── download_genomes.py      # 基因组下载 / Genome downloading
└── data/                        # 数据目录 / Data directory
    └── assembly_summary.txt     # 基因组信息文件 / Genome information file
```

## 🚀 快速开始 / Quick Start

### 1. 环境准备 / Environment Setup

```bash
# 安装Python依赖 / Install Python dependencies
pip install pandas biopython numpy

# 安装外部工具 / Install external tools
conda install -c bioconda prokka roary mafft primer3
```

### 2. 基本使用 / Basic Usage

```bash
# 运行完整流水线 / Run complete pipeline
python run_primer_design.py --genus Terrisporobacter --threads 8

# 从特定步骤开始 / Start from specific step
python run_primer_design.py --genus Terrisporobacter --start-step 3

# 设置最大基因数量 / Set maximum genes
python run_primer_design.py --genus Terrisporobacter --max-genes 20
```

### 3. 独立工具使用 / Independent Tool Usage

```bash
# 引物质量评分 / Primer quality ranking
python utils/rank_primers.py input_primers.csv output_ranked.csv --threads 4

# Primer3结果解析 / Primer3 result parsing
python utils/parse_primer3.py primer3_output.txt parsed_primers.csv

# 基因组下载 / Genome downloading
python utils/download_genomes.py --genus Terrisporobacter --output genomes/
```

## 🔧 配置管理 / Configuration Management

项目使用配置文件 `.project_config` 管理全局编码规范：

The project uses `.project_config` file to manage global coding standards:

```python
from primer_design_toolkit.config_manager import get_config_manager

# 获取配置管理器 / Get configuration manager
config = get_config_manager()

# 检查是否使用双语注释 / Check if bilingual comments should be used
if config.should_use_bilingual_comments():
    # 使用双语注释格式 / Use bilingual comment format
    comment_format = config.get_comment_format()

# 检查是否使用英文输出 / Check if English output should be used
if config.should_use_english_output():
    print("Using English output")
```

## 📊 流水线步骤 / Pipeline Steps

### 步骤1: 基因组下载 / Step 1: Genome Download
- 基于assembly_summary.txt下载目标和外群基因组
- Download target and outgroup genomes based on assembly_summary.txt

### 步骤2: Prokka注释 / Step 2: Prokka Annotation  
- 并行进行基因组注释，生成GFF文件
- Parallel genome annotation to generate GFF files

### 步骤3: Roary泛基因组分析 / Step 3: Roary Pan-genome Analysis
- 构建泛基因组，识别核心和附属基因
- Build pan-genome and identify core and accessory genes

### 步骤4: 特异性基因筛选 / Step 4: Specific Gene Filtering
- 筛选只在目标菌属中存在的基因
- Filter genes present only in target genus

### 步骤5: 基因序列提取 / Step 5: Gene Sequence Extraction
- 提取特异性基因的DNA序列
- Extract DNA sequences of specific genes

### 步骤6: 多序列比对 / Step 6: Multiple Sequence Alignment
- 使用MAFFT进行多序列比对，识别保守区域
- Use MAFFT for multiple sequence alignment and identify conserved regions

### 步骤7: 引物设计 / Step 7: Primer Design
- 使用Primer3设计引物对
- Use Primer3 to design primer pairs

### 步骤8: 结果解析 / Step 8: Result Parsing
- 解析Primer3输出，提取引物信息
- Parse Primer3 output and extract primer information

### 步骤9: 质量评估 / Step 9: Quality Assessment
- 基于10维评分体系评估引物质量
- Assess primer quality based on 10-dimensional scoring system

## 🎯 质量评分体系 / Quality Scoring System

引物质量评分基于以下10个维度：

Primer quality scoring based on 10 dimensions:

1. **产物大小 Product Size** (权重 Weight: 10%)
2. **熔解温度 Melting Temperature** (权重 Weight: 15%)  
3. **Tm差异 Tm Difference** (权重 Weight: 10%)
4. **GC含量 GC Content** (权重 Weight: 15%)
5. **GC差异 GC Difference** (权重 Weight: 5%)
6. **引物长度 Primer Length** (权重 Weight: 10%)
7. **自身二聚体 Self Dimer** (权重 Weight: 10%)
8. **末端二聚体 End Dimer** (权重 Weight: 10%)
9. **发夹结构 Hairpin Structure** (权重 Weight: 5%)
10. **引物对互补性 Primer Pair Complementarity** (权重 Weight: 10%)

### 质量等级 / Quality Grades
- **A+ (90-100分)**: 优秀引物 / Excellent primers
- **A (85-89分)**: 良好引物 / Good primers  
- **B+ (80-84分)**: 可用引物 / Usable primers
- **B (75-79分)**: 一般引物 / Fair primers
- **C+ (70-74分)**: 较差引物 / Poor primers
- **C (65-69分)**: 差引物 / Bad primers
- **D (<65分)**: 不推荐 / Not recommended

## 🔍 代码质量检查 / Code Quality Check

项目提供自动化代码质量检查：

The project provides automated code quality checking:

```python
from primer_design_toolkit.config_manager import ProjectConfigManager

config = ProjectConfigManager()

# 验证代码标准 / Validate code standards
code_text = """
def example_function():
    # 示例函数 / Example function
    print("Hello World!")
"""

results = config.validate_code_standards(code_text)
print(f"Validation results: {results}")
```


## 📄 许可证 / License

MIT License - 详见LICENSE文件 / See LICENSE file for details


>>>>>>> 45c61a66b7f92853ff8769c02cbf6476122f7d17
