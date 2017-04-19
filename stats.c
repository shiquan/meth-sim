#include "utils.h"
//#include "stats.h"
#include <time.h>
#include <math.h>

static int seed = -1;

void update_seed(int s)
{
    seed = s;
}

int get_seed()
{
    if ( seed <= 0 )
        return time(0)&0x7fffffff;
    return seed;
}

// Box muller method
double ran_normal(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

// random return 1 or 0 based on the probably 
int random_box1(float prob)
{
    return (rand()/(double)RAND_MAX) < prob;    
}

int random_box2(float prob)
{
    float r = ran_normal(0.5, 0.25);
    while ( r < 0 )
        r = 1 + r;
    while  ( r > 1)
        r = r - 1;
    return r < prob;
}


#ifdef STAT_MAIN

int main()
{
    float f;
    int i;
    float total = 0;
    for ( f = 0.001; f <= 1.0; f+=0.001) {
        int c1;
        int c2;
        float p1, p2;
        for (c1=0, c2=0, i = 0; i < 100; ++i ) {
            int t1 = random_box1(f);
            if ( t1 ) {
                c1++;
            } else {
                c2++;
            }            
        }
        p1 = (float)c1/(c1+c2);
        for (c1=0, c2=0, i = 0; i < 100; ++i ) {
            int t1 = random_box2(f);
            if ( t1 ) {
                c1++;
            } else {
                c2++;
            }            
        }
        p2 = (float)c1/(c1+c2);
        printf("%f\t%f\t%f\t%d\n", f, p1, p2, (int)(-1000*log(f)));
    }
    return 0;
}
#endif
