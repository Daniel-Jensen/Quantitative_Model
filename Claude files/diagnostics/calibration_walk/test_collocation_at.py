# RUN test_collocation's real solve at a chosen pd_rest, to tell a pre-existing failure from one the 0.02% change made.
# usage: python3 test_collocation_at.py 0.001
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path[:0] = [ROOT, os.path.join(ROOT, "Claude files", "global_projection", "tests")]

if __name__ == "__main__":
    import test_collocation as T
    get = T.get_ss

    def patched():
        cal, ss = get()
        cal["pd_rest"] = float(sys.argv[1])
        return cal, ss
    T.get_ss = patched
    T.test_solve_reaches_the_floor()
    print("PASSED at pd_rest", sys.argv[1])
