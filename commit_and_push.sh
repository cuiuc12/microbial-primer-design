#!/bin/bash

# 🔄 提交精简版本到GitHub
echo "🔄 Committing cleaned up version to GitHub..."

# 检查是否在正确目录
if [[ ! -f "run_primer_design.py" ]]; then
    echo "❌ Error: Please run this script from the project root directory"
    exit 1
fi

# 检查Git状态
echo "📋 Current Git status:"
git status

echo ""
read -p "📝 Continue with commit? (y/n): " -n 1 -r
echo ""

if [[ $REPLY =~ ^[Yy]$ ]]; then
    # 添加所有更改
    echo "📥 Adding all changes..."
    git add -A
    
    # 提交更改
    echo "💾 Committing changes..."
    git commit -m "🧹 Clean up project: Remove docs and config files, keep only core functionality

- Remove unnecessary documentation files (CODING_STANDARDS.md, INSTALLATION_GUIDE.md, etc.)
- Remove config manager and project config files  
- Simplify dependency files (requirements.txt, environment.yml)
- Streamline README.md with only core usage instructions
- Project now contains only essential files for running

Core retained files:
- run_primer_design.py (main script)
- primer_design_toolkit/ (core toolkit)
- utils/ (standalone tools)
- Basic config files (requirements.txt, .gitignore, LICENSE)

Clean, focused, and ready to use!"
    
    # 推送到GitHub
    echo "🚀 Pushing to GitHub..."
    git push origin main
    
    if [[ $? -eq 0 ]]; then
        echo "✅ Successfully pushed to GitHub!"
        echo "🎯 Your repository is now clean and streamlined!"
    else
        echo "❌ Push failed. Please check your authentication and try again."
        echo "💡 You may need to use your Personal Access Token instead of password."
    fi
else
    echo "❌ Operation cancelled."
fi 