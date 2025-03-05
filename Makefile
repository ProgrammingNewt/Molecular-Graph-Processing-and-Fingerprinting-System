install:
	python -m pip install networkx pytest matplotlib

test:
	pytest molTesting.py

clean:
	rm -rf __pycache__ .pytest_cache