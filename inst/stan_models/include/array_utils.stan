//functions {
    // count distinct elements of ixs
    int ndistinct(int[] ixs, int maxix) {
      int ix2used[maxix] = rep_array(0, maxix);
      for (ix in ixs) {
        ix2used[ix] = 1;
      }
      return sum(ix2used);
    }

    // distinct elements of ixs
    int[] distinct(int[] ixs, int maxix) {
      int ixs_pos[maxix] = rep_array(0, maxix);
      int ixs_distinct[ndistinct(ixs, maxix)] = rep_array(0, maxix);
      int offset = 0;
      for (ix in ixs) {
        if (ixs_pos[ix] == 0) {
            offset += 1;
            ixs_distinct[offset] = ix;
            ixs_pos[ix] = offset;
        }
      }
      return ixs_distinct;
    }

    // count how many type each group ix occurred
    int[] group_sizes(int[] groupixs, int ngroups) {
        int counts[ngroups] = rep_array(0, ngroups);
        for (i in 1:num_elements(groupixs)) {
          counts[groupixs[i]] += 1;
        }
        return counts;
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