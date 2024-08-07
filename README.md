# simplex

A [Simplex algorithm](https://en.wikipedia.org/wiki/Simplex_algorithm) calculator that returns certificates alongside the outome of an [linear program](https://en.wikipedia.org/wiki/Linear_programming) in standard equality form.


## How to use

Clone the repo and navigate to its root.

Now, run the following commands:
```
pip install -r requirements.txt
```

Update the values of `A`, `b`, and `c` to reflect your linear program. It should be of the form:

$$\text{max } c^T x\\ \text{s.t. } Ax = b\\ x \ge 0$$


Then just run:
```
python simplex.py
```
