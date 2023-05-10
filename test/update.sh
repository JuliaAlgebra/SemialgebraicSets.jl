find . -type f -name "*.txt" -print0 | xargs -0 sed -i '' -e 's/foo/bar/g'

