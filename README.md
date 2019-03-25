# NeedlemanWunsh-SmithWaterman

Implementation of NeedlemanWunsh-SmithWaterman algorithm in python.
The file align.py takes in two arguments on the command line, the first being the input file, the second being the output filename.

Example input is located in "examples/alignment_example.input.annotated.txt"
Example output is located in "examples/alignment_example.output.txt"

Output is ranked in order of best score, with the best score indicated at the top of the file.

Test data is located in the folder "test_data", along with the results from when the algorithm is run on the corresponding input.

To run the script, input the following command:
```python
python align.py alignment0.input.txt alignment0.output.txt
```

