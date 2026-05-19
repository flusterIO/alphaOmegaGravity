copy_notebook:
	cp -r ~/Desktop/notes/content/physics/pseudo_notebooks/python_package_as_notebook/ {{justfile_directory()}}/physics/notebook

release: copy_notebook
	python -m build && python -m twine upload dist/*
