#!/bin/bash

# 🔄 重新初始化Git仓库并提交到GitHub
echo "🔄 Resetting Git repository and committing fresh to GitHub..."

# 检查是否在正确目录
if [[ ! -f "run_primer_design.py" ]]; then
    echo "❌ Error: Please run this script from the project root directory"
    exit 1
fi

echo "⚠️  Warning: This will delete all Git history and create a fresh repository!"
read -p "📝 Continue? (y/n): " -n 1 -r
echo ""

if [[ $REPLY =~ ^[Yy]$ ]]; then
    # 删除现有的.git目录
    echo "🗑️ Removing existing Git history..."
    rm -rf .git
    
    # 重新初始化Git仓库
    echo "🆕 Initializing new Git repository..."
    git init
    git branch -M main
    
    # 添加所有文件
    echo "📥 Adding all files..."
    git add .
    
    # 创建初始提交
    echo "💾 Creating initial commit..."
    git commit -m "🎯 Microbial Primer Design Toolkit - Clean Release

A comprehensive pipeline toolkit for microbial-specific primer design.

Features:
✅ Complete 9-step automated pipeline
✅ Parallel processing for 89% of workflows  
✅ 3-4x performance improvement
✅ 10-dimensional primer quality scoring
✅ Support for NCBI genome data
✅ Easy one-command usage

Core Components:
- run_primer_design.py (main script)
- primer_design_toolkit/ (core modules)
- utils/ (standalone tools)
- Complete documentation and examples

Ready to use for microbial-specific primer design!"
    
    # 添加远程仓库
    echo "🔗 Adding remote repository..."
    git remote add origin https://github.com/cuiuc12/microbial-primer-design.git
    
    # 推送到GitHub (强制)
    echo "🚀 Pushing to GitHub..."
    git push -u origin main --force
    
    if [[ $? -eq 0 ]]; then
        echo ""
        echo "✅ Successfully pushed to GitHub!"
        echo "🎯 Your repository is now clean and ready!"
        echo "📍 Repository URL: https://github.com/cuiuc12/microbial-primer-design"
        echo ""
        echo "🔗 You can now share this repository with others!"
    else
        echo "❌ Push failed. Please check your authentication."
        echo "💡 Make sure you're using your Personal Access Token."
    fi
else
    echo "❌ Operation cancelled."
fi 