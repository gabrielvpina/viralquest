# 1. Create a branch to preserve the old main
git checkout main
git checkout -b legacy-v2

# 2. Push it to remote so it's preserved there too
git push origin legacy-v2

# 3. Reset main to total-refactor
git checkout main
git reset --hard total-refactor

# 4. Force push the new main
git push origin main --force
