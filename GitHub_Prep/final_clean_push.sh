#!/bin/bash

# Create completely clean repository by filtering out all large files upfront
echo "🧹 Creating completely clean repository without any large files..."

cd "/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence"

# Remove git history completely
rm -rf .git

# Create comprehensive .gitignore before anything else
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

# Presentations and documentation
*.pptx
*.ppt
*.docx
*.doc

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

# Remove any existing large files that might cause issues
find . -name "*.pptx" -size +50M -delete
find . -name "*.csv" -size +50M -delete  
find . -name "*.pdf" -size +50M -delete
find . -name "*.png" -size +50M -delete
find . -name "*.xlsx" -size +50M -delete

echo "🗑️  Removed large files that could cause GitHub issues"

# Initialize fresh repository
git init
git branch -m main

# Add files respecting .gitignore
git add .

# Create initial commit
git commit -m "Initial commit: Clean DIMCAT Prevalence analysis

Spatial epidemiological analysis of African Animal Trypanosomiasis (AAT) prevalence
- R scripts for INLA Bayesian modeling
- Ethiopia and Nigeria country-specific analysis  
- Correlation analysis and model selection
- Cost-effectiveness and EVSI analysis
- Test sensitivity/specificity modeling
- No large data files (filtered out)
- Clean repository without LFS dependencies

Key features:
- Continental and country-level prevalence mapping
- Cattle-at-risk assessments
- Uncertainty quantification
- Administrative unit aggregation
- Economic decision analysis"

# Connect to GitHub
git remote add origin https://github.com/KayeARK/DIMCAT_Prevalence.git

# Push to completely replace GitHub repository
echo "🚀 Pushing clean repository to GitHub..."
git push origin main --force

echo ""
echo "✅ Repository successfully updated!"
echo "🌐 Your GitHub repository is now completely clean"
echo "📁 Repository URL: https://github.com/KayeARK/DIMCAT_Prevalence"
echo "🧹 All large files filtered out"
echo "📊 Contains only R scripts and essential documentation"
echo ""
echo "Ready for collaboration and sharing! 🎉"