//functions {
    int ndistinct(int[] ixs, int maxix) {
      int ix2used[maxix] = rep_array(0, maxix);
      for (ix in ixs) {
        ix2used[ix] = 1;
      }
      return sum(ix2used);
    }

    int[] one_to(int n) {
      int res[n];
      for (i in 1:n) res[i] = i;
      return res;
    }

    // return the ordering of vals so that nonzeros are put first
    int[] nonzeros_first(int[] vals) {
        int two_or_one[size(vals)];
        for (i in 1:size(vals)) {
            two_or_one[i] = vals[i] > 0 ? 1 : 2;
        }
        return sort_indices_asc(two_or_one);
    }
//}