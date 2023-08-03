#include <math.h>
#include <stdint.h>


struct acrGrid {
        double bin_width;
};

/*                                 \                                          */
/*                                  \                                    __-- */
/*                                   \                               __--     */
/*   +---------+---------+---------+--x------+---------+---------+o--------+  */
/*   |         |         |         |   \     |         |    __-- |         |  */
/*   |         |         |         |    \    |         |__--     |         |  */
/*   |         |         |         |     \   |      __-|         |         |  */
/*   |         |         |         |      \  |  __--   |         |         |  */
/*   +----|----+----|----+----|----+----|--x-+o---|----+----|----+----|----+  */
/*                                                                            */
/*       -3        -2        -1         0         1         2         3       */


int64_t acr_query_bin(double x, double bin_width)
{
        double bin = x / bin_width;
        return (int64_t)round(bin);
}

int acr_query_bins(
        double bin_width,
        double bin_height,
        double ray_x,
        double ray_y,
        double ray_cx,
        double ray_cy,
        uint64_t intersections_capacity,
        int64_t *intersections_x,
        int64_t *intersections_y,
        uint64_t *intersections_num)
{
        (*intersections_num) = 0u;



        return 0;
}

