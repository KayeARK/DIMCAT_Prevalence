#!/bin/bash

# Reset existing GitHub repository with clean history (no LFS)
# This completely replaces your GitHub repository content while keeping the same URL

echo "🔄 Updating your existing GitHub repository with clean history..."
echo "⚠️  This will replace the GitHub repository history but keep the same URL"
echo ""

cd "/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence"

# Create a backup of current directory
echo "📦 Creating backup..."
cp -r . ../DIMCAT_Prevalence_backup_$(date +%Y%m%d_%H%M%S)

# Create completely new git repository
echo "🗑️  Removing all git history..."
rm -rf .git

echo "🆕 Creating fresh git repository..."
git init
git branch -m main

# Make sure we have a comprehensive .gitignore
echo "📝 Updating .gitignore..."
cat > .gitignore << 'EOF'
# Large data files
*.csv
*.xlsx
*.xls
*.pdf
*.png
*.jpg
*.jpeg
*.gif
*.tif
*.tiff
*.rds
*.RData
*.Rdata

# Results and outputs
Results/
results/
output/
outputs/
plots/
figures/

# R specific
.Rproj.user/
.Rhistory
.RData
.Ruserdata
*.Rproj

# System files
.DS_Store
.DS_Store?
._*
.Spotlight-V100
.Trashes
ehthumbs.db
Thumbs.db

# IDE
.vscode/
.idea/

# Temporary files
*~
*.tmp
*.temp

# Log files
*.log

# Archives
*.zip
*.tar.gz
*.rar
EOF

# Add all files respecting .gitignore
echo "📁 Adding files to repository..."
git add .

# Create initial commit
echo "💾 Creating initial commit..."
git commit -m "Initial commit: Clean DIMCAT Prevalence analysis codebase

- R scripts for spatial epidemiological analysis
- INLA Bayesian modeling for African Animal Trypanosomiasis
- Ethiopia and Nigeria country-specific analysis
- Correlation analysis and model selection
- No large data files (see .gitignore)
- Clean repository without LFS dependencies"

# Add existing GitHub remote
echo "🔗 Connecting to existing GitHub repository..."
git remote add origin https://github.com/KayeARK/DIMCAT_Prevalence.git

# Force push to completely replace GitHub repository
echo "🚀 Updating GitHub repository..."
echo "⚠️  This will completely replace the existing repository content on GitHub"
git push origin main --force

echo ""
echo "✅ Repository updated successfully!"
echo "🌐 Your GitHub repository has been completely refreshed"
echo "📂 Backup created at: ../DIMCAT_Prevalence_backup_$(date +%Y%m%d_%H%M%S)"
echo "🔗 Repository URL remains: https://github.com/KayeARK/DIMCAT_Prevalence"
echo ""
echo "The repository now contains:"
echo "- All your R scripts and code"
echo "- Clean git history (no LFS issues)"
echo "- Comprehensive .gitignore for large files"
echo "- Ready for collaboration and sharing"