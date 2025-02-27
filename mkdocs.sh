cd docs && poetry run sphinx-build -b html -E source ../gh-pages
cd ..
echo open gh-pages/index.html