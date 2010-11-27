/*
  Stats.java
  Copyright © 2010 David M. Anderson

  Statitistical estimates and hypothesis tests.
  NOTES:
  1. Variance() returns the unbiased estimate, dividing by (N-1), rather than
     by N, as is sometimes done.
  2. Median() modifies the sample data only to the extent that it reorders
     the elements.
  3. Hypothesis tests return the probability of the observed result, given
     the (null) hypothesis. So a smaller value suggests rejection.
  4. MeanTest() tests the hypothesis that the population mean is equal to the
     hypothetical mean. It is valid for a normal distribution or a large
     sample size. Likewise for VarianceTest().
  5. KolmogorovSmirnovTest() tests the hypothesis that the population
     distribution function is that given by the hypothDist. It sorts the
     sample data if not already sorted.
  6. ChiSquareGoodnessOfFitTest() tests whether binned or categorized data
     comes from a distribution with the hypothesized frequencies. If, as is
     commonly the case, the hypothesis specifies probabilities, leave
     probabilities = true. If, however, the expected frequencies are directly
     known, set probabilities=false. If other parameters of the expected
     distribution were determined from the data, constraints should be set
     to the number of such free parameters have been fitted.
  7. MeansTest() tests the hypothesis that the means of two populations are
     the same. It is valid for normal distributions or large sample sizes.
     Set equalVariances = true if the variances of both distributions are
     believed to be the same; false if they are significantly different.
     VariancesTest() can be used to decide this.
  8. MediansTest() tests the hypothesis that the medians of two populations
     are equal.
  9. The second form of KolmogorovSmirnovTest() tests the hypothesis that
     two samples are drawn from populations with the same distribution.
  10. ContingencyTableTest() 
*/

package us.EpsilonDelta.math;

import java.util.*;


//*****************************************************************************


public class Stats
{                                                                       //Stats
//-----------------------------------------------------------------------------

    public static
    double 
    mean( final List< Double > sample )
    {
        if ( sample.size() <= 0 )
            throw new IllegalArgumentException( "mean: sample size <= 0" );
        double sum = 0.;
        for ( Double d : sample )
            sum += d;
        return sum / sample.size();
    }

//-----------------------------------------------------------------------------

    public static
    double 
    variance( final List< Double > sample, double mean )
    {
        int sampleSize = sample.size();
        if ( sampleSize <= 1 )
            throw new IllegalArgumentException( "variance: sample size <= 1" );
        //Chan, T.F., et al., American Statistician, vol 37 (1983), p 242-47.
        double diffSum = 0.;
        double diffSqrSum = 0.;
        for ( Double d : sample )
        {
            double diff = d - mean;
            diffSum += diff;
            diffSqrSum += diff * diff;
        }
        return ( diffSqrSum  -  (diffSum * diffSum) / sampleSize )
                / (sampleSize - 1.);
    }

//-----------------------------------------------------------------------------

    public static
    double 
    median( List< Double > sample )
    {
        return median( sample, false );
    }
    
//.............................................................................
    
    public static
    double 
    median( List< Double > sample, boolean alreadySorted )
    {
        int sampleSize = sample.size();
        if ( sampleSize <= 0 )
            throw new IllegalArgumentException( "median: sample size <= 0" );
        int n2 = (sampleSize - 1) / 2;
        if ( ! alreadySorted )
            Collections.sort( sample );
        double median = sample.get( n2 );
        if ( (sampleSize & 1) == 0 ) //even
        {
            median += sample.get( n2 + 1 );
            median /= 2.;
        }
        return median;
    }

//=============================================================================

    public static
    TTestResult
    meanTest( int sampleSize, double sampleMean, double sampleVariance,
              double hypothMean, int tails )
    {
        if ( sampleSize <= 1 )
            throw new IllegalArgumentException( "meanTest: sampleSize <= 1" );
        if ( sampleVariance < 0. )
            throw new IllegalArgumentException( "meanTest: "
                                                + "sampleVariance < 0" );
        double prob = 0.;
        double t = 0.;
        int degreesOfFreedom = sampleSize - 1;
        if ( sampleVariance == 0. )
        {
            prob = (sampleMean == hypothMean) ? 1. : 0.;
            t = (sampleMean == hypothMean) ? 0.
                    : (sampleMean < hypothMean) ? Double.NEGATIVE_INFINITY
                    : Double.POSITIVE_INFINITY;
        }
        else
        {
            t = (sampleMean - hypothMean) * Math.sqrt( (double) sampleSize )
                    / Math.sqrt( sampleVariance );
            prob = ProbDist.t_DF( t, degreesOfFreedom );
            if ( t > 0 )
                prob = 1. - prob;
            if ( tails == 2 )
                prob *= 2.;
        }
        return new TTestResult( prob, t, degreesOfFreedom );
    }

//-----------------------------------------------------------------------------

    public static
    ChiSquareTestResult
    varianceTest( int sampleSize, double sampleVariance, double hypothVariance,
                  int tails )
    {
        if ( sampleSize <= 1 )
            throw new IllegalArgumentException( "varianceTest: "
                                                + "sampleSize <= 1" );
        if ( sampleVariance < 0. )
            throw new IllegalArgumentException( "varianceTest: "
                                                + "sampleVariance < 0" );
        if ( hypothVariance < 0. )
            throw new IllegalArgumentException( "varianceTest: "
                                                + "hypothVariance < 0" );
        double prob = 0.;
        double chiSquare = 0.;
        int degreesOfFreedom = sampleSize - 1;
        if ( hypothVariance == 0. )
        {
            prob = (sampleVariance == 0.) ? 1. : 0.;
            chiSquare = (sampleVariance == 0.) ? 0. : Double.POSITIVE_INFINITY;
        }
        else
        {
            chiSquare = sampleVariance * (sampleSize - 1.) / hypothVariance;
            prob = ProbDist.chiSquare_DF( chiSquare, degreesOfFreedom );
        }
        if ( sampleVariance > hypothVariance )
            prob = 1. - prob;
        if ( tails == 2 )
            prob *= 2.;
        return new ChiSquareTestResult( prob, chiSquare, degreesOfFreedom );
    }

//=============================================================================

    public static
    ChiSquareTestResult
    chiSquareGoodnessOfFitTest( final List< Integer > sampleFreqs,
                                final List< Double > hypothFreqs,
                                boolean probabilities, int constraints )
    {
        int numBins = sampleFreqs.size();
        if ( hypothFreqs.size() != numBins )
            throw new IllegalArgumentException( "chiSquareGoodnessOfFitTest"
                                                + "sampleFreqs size "
                                                + "!= hypothFreqs size" );
        int degreesOfFreedom = numBins - constraints;
        int sampleSize = 0;
        if ( probabilities )
        {
            --degreesOfFreedom;
            for ( int f : sampleFreqs )
                sampleSize += f;
        }
        double chiSquare = 0.;
        for ( int i = 0; i < numBins; ++i )
        {
            double hypFreq = hypothFreqs.get( i );
            double sampleFreq = sampleFreqs.get( i );
            if ( probabilities )
                hypFreq *= sampleSize;
            if ( (hypFreq == 0.) && (sampleFreq > 0) )
            {
                return new ChiSquareTestResult( 0., Double.POSITIVE_INFINITY,
                                                degreesOfFreedom );
            }
            if ( (hypFreq < 0.5) && (sampleFreq == 0) )
            {
                --degreesOfFreedom;
                continue;
            }
            double diff = sampleFreq - hypFreq;
            chiSquare += diff * diff / hypFreq;
        }
        double prob = 1. - ProbDist.chiSquare_DF( chiSquare, degreesOfFreedom );
        return new ChiSquareTestResult( prob, chiSquare, degreesOfFreedom );
    }

//=============================================================================

    public static
    TTestResult
    meansTest( int sampleSize1, double sampleMean1, double sampleVariance1,
               int sampleSize2, double sampleMean2, double sampleVariance2,
               boolean equalVariances, int tails )
    {
        if ( sampleSize1 <= 1 )
            throw new IllegalArgumentException( "meansTest: sampleSize1 <= 1" );
        if ( sampleVariance1 < 0. )
            throw new IllegalArgumentException( "meansTest: "
                                                + "sampleVariance1 < 0" );
        if ( sampleSize2 <= 1 )
            throw new IllegalArgumentException( "meansTest: sampleSize2 <= 1" );
        if ( sampleVariance2 < 0. )
            throw new IllegalArgumentException( "meansTest: "
                                                + "sampleVariance2 < 0" );
        double pooledVariance;
        int degreesOfFreedom;
        if ( equalVariances )
        {
            pooledVariance = ((sampleVariance1 * (sampleSize1 - 1.)
                               +  sampleVariance2 * (sampleSize2 - 1.))
                              / (sampleSize1 + sampleSize2 - 2.))
                    *  ((1. / sampleSize1) + (1. / sampleSize2));
            degreesOfFreedom = sampleSize1 + sampleSize2 - 2;
        }
        else
        {
            double vn1 = sampleVariance1 / sampleSize1;
            double vn2 = sampleVariance2 / sampleSize2;
            pooledVariance = vn1 + vn2;
            double dof = ( ((vn1 + vn2) * (vn1 + vn2))
                           / ((vn1 * vn1 / (sampleSize1 - 1.))
                              +  (vn2 * vn2 / (sampleSize2 - 1.))) );
            degreesOfFreedom = (int) dof;
        }
        double prob = 0.;
        double t = 0.;
        assert pooledVariance >= 0.;
        if ( pooledVariance == 0. )
        {
            prob = (sampleMean1 == sampleMean2) ? 1. : 0.;
            t = (sampleMean1 == sampleMean2) ? 0.
                    : (sampleMean1 < sampleMean2) ? Double.NEGATIVE_INFINITY
                    : Double.POSITIVE_INFINITY;
        }
        else
        {
            t = (sampleMean1 - sampleMean2) / Math.sqrt( pooledVariance );
            prob = ProbDist.t_DF( t, degreesOfFreedom );
            if ( t > 0 )
                prob = 1. - prob;
            if ( tails == 2 )
                prob *= 2.;
        }
        return new TTestResult( prob, t, degreesOfFreedom );
    }

//-----------------------------------------------------------------------------

    public static
    FTestResult
    variancesTest( int sampleSize1, double sampleVariance1,
                   int sampleSize2, double sampleVariance2,
                   int tails )
    {
        if ( sampleSize1 <= 1 )
            throw new IllegalArgumentException( "variancesTest: "
                                                + "sampleSize1 <= 1" );
        if ( sampleVariance1 < 0. )
            throw new IllegalArgumentException( "variancesTest: "
                                                + "sampleVariance1 < 0" );
        if ( sampleSize2 <= 1 )
            throw new IllegalArgumentException( "variancesTest: "
                                                + "sampleSize2 <= 1" );
        if ( sampleVariance2 < 0. )
            throw new IllegalArgumentException( "variancesTest: "
                                                + "sampleVariance2 < 0" );
        double prob = 0.;
        double f = 0.;
        int dof1 = (sampleSize1 - 1);
        int dof2 = (sampleSize2 - 1);
        if ( sampleVariance2 == 0. )
        {
            prob = (sampleVariance1 == 0.) ? 1. : 0.;
            f = (sampleVariance1 == 0.) ? 0. : Double.POSITIVE_INFINITY;
        }
        else
        {
            f = sampleVariance1 / sampleVariance2;
            prob = ProbDist.F_DF( f, dof1, dof2 );
            if ( f > 1. )
                prob = 1. - prob;
            if ( tails == 2 )
                prob *= 2.;
        }
        return new FTestResult( prob, f, dof1, dof2 );
    }

//=============================================================================

    public static
    double 
    mediansTest( final List< Double > sample1,
                 final List< Double > sample2, int tails )
    {
        List< Double > grandSample = new ArrayList< Double >();
        grandSample.addAll( sample1 );
        grandSample.addAll( sample2 );
        double grandMedian = median( grandSample );
        int numBelowMedian = 0;
        int numAboveMedian = 0;
        for ( Double d : sample1 )
        {
            if ( d < grandMedian )
                ++numBelowMedian;
            else if ( d > grandMedian )
                ++numAboveMedian;
        }
        int grandSize = grandSample.size();
        int sampleSize1 = sample1.size();
        int x = Math.min( numBelowMedian, numAboveMedian );
        double prob = ProbDist.hypergeometric_DF( x, grandSize, grandSize / 2,
                                                  sampleSize1 );
        if ( tails == 2 )
            prob *= 2.;
        return prob;
    }

//=============================================================================

    public static
    double 
    kolmogorovSmirnovTest( List< Double > sample1, List< Double > sample2,
                           boolean alreadySorted )
    {
        if ( ! alreadySorted )
        {
            Collections.sort( sample1 );
            Collections.sort( sample2 );
        }
        double d = 0.;
        int sampleSize1 = sample1.size();
        int sampleSize2 = sample2.size();
        if ( sampleSize1 == 0 )
            throw new IllegalArgumentException(
                "kolmogorovSmirnovTest: sampleSize1 = 0" );
        if ( sampleSize2 == 0 )
            throw new IllegalArgumentException(
                "kolmogorovSmirnovTest: sampleSize2 = 0" );
        double recipSize1 = 1. / sampleSize1;
        double recipSize2 = 1. / sampleSize2;
        double DF1 = 0.;
        double DF2 = 0.;
        int i1 = 0;
        int i2 = 0;
        while ( (i1 < sampleSize1) && (i2 < sampleSize2) )
        {
            double s1 = sample1.get( i1 );
            double s2 = sample2.get( i2 );
            if ( s1 <= s2 )
            {
                ++i1;
                DF1 += recipSize1;
            }
            if ( s2 <= s1 )
            {
                ++i2;
                DF2 += recipSize2;
            }
            double diff = Math.abs( DF1 - DF2 );
            if ( diff > d )
                d = diff;
        }
        double sqrtSize = Math.sqrt( (double) sampleSize1 * sampleSize2 )
                / (sampleSize1 + sampleSize2);
        double x = d * (sqrtSize  +  0.12  +  0.11 / sqrtSize);
        return 1. - ProbDist.kolmogorovSmirnov_DF( x );
    }

//=============================================================================

    public static
    CorrelationTestResult
    linearCorrelationTest( final List< double[] > samples, int tails )
    {
        int sampleSize = samples.size();
        if ( sampleSize < 3 )
            throw new IllegalArgumentException( "linearCorrelationTest: "
                                                + "sampleSize < 3" );
        double mean0 = 0.;
        double mean1 = 0.;
        for ( double[] da : samples )
        {
            if ( da.length != 2 )
                throw new IllegalArgumentException( "linearCorrelationTest: "
                                                    + "samples[].size != 2" );
            mean0 += da[0];
            mean1 += da[1];
        }
        mean0 /= sampleSize;
        mean1 /= sampleSize;
        double sum0 = 0.;
        double sum1 = 0.;
        double sum00 = 0.;
        double sum11 = 0.;
        double sum01 = 0.;
        for ( double[] da : samples )
        {
            double diff0 = da[0] - mean0;
            double diff1 = da[1] - mean1;
            sum0 += diff0;
            sum1 += diff1;
            sum00 += diff0 * diff0;
            sum11 += diff1 * diff1;
            sum01 += diff0 * diff1;
        }
        //Chan, T.F., et al., American Statistician, vol 37 (1983), p 242-47.
        double variance0 = ( sum00  -  (sum0 * sum0) / sampleSize )
                / (sampleSize - 1.);
        double variance1 = ( sum11  -  (sum1 * sum1) / sampleSize )
                / (sampleSize - 1.);
        int degreesOfFreedom = sampleSize - 2;
        double covariance = ( sum01 - (sum0 * sum1) / sampleSize )
                / (sampleSize - 1.);
        double denom = Math.sqrt( variance0 * variance1 );
        double r = 0.;
        double t = 0.;
        double prob = 0.;
        if ( denom == 0. )
            r = (covariance < 0) ? -1. : 1.;
        else
            r = covariance / denom;
        if ( r <= -1. )
        {
            t = Double.NEGATIVE_INFINITY;
            prob = 0.;
        }
        else if ( r >= 1. )
        {
            t = Double.POSITIVE_INFINITY;
            prob = 0.;
        }
        else
        {
            t = r * Math.sqrt( (sampleSize - 2.) / (1. - r * r) );
            prob = ProbDist.t_DF( t, degreesOfFreedom );
            if ( t > 0. )
                prob = 1. - prob;
            if ( tails == 2 )
                prob *= 2.;
        }
        return new CorrelationTestResult(
            new TTestResult( prob, t, degreesOfFreedom ),
            r, mean0, mean1, variance0, variance1 );
    }

//=============================================================================

    public static
    SimpleRegressionResult
    simpleLinearRegression( int sampleSize,
                            double meanX, double meanY,
                            double varianceX, double varianceY,
                            double pearsonsR )
    {
        if ( sampleSize <= 2 )
            throw new IllegalArgumentException( "simpleLinearRegression: "
                                                + "sampleSize <= 2" );
        if ( varianceX <= 0. )
            throw new IllegalArgumentException( "simpleLinearRegression: "
                                                + "varianceX <= 0" );
        if ( varianceY < 0. )
            throw new IllegalArgumentException( "simpleLinearRegression: "
                                                + "varianceY < 0" );
        double beta = pearsonsR * Math.sqrt( varianceY / varianceX );
        double alpha = meanY - beta * meanX;
        double residualVariance = ((sampleSize - 1.) / (sampleSize - 2.))
                * varianceY * (1. - pearsonsR * pearsonsR);
        double varianceAlpha = residualVariance
                * ((1. / sampleSize)
                   +  (meanX * meanX) / ((sampleSize - 1.) * varianceX));
        double varianceBeta = residualVariance / ((sampleSize - 1) * varianceX);
        return new SimpleRegressionResult( alpha, beta, varianceAlpha,
                                           varianceBeta, residualVariance );
    }
    
//=============================================================================

/*!!!
  namespace
  {                                                                   //namespace

  class CompareArrays2
  {
  public:
  CompareArrays2( int index )
  :   m_index( index )
  {
  }

  boolean operator()( final array< double, 2 > & arr1,
  final array< double, 2 > & arr2 )
  {
  return ( arr1[ m_index ] < arr2[ m_index ] );
  }

  private:
  int m_index;
  };

  }                                                                   //namespace

*/

//.............................................................................

    public static
    double 
    spearmansRankCorrelationTest( final List< double[] > samples,
                                  Double spearmansR, Double t )
    {
/*!!!
  int sampleSize = samples.size();
  if ( sampleSize <= 3 )
  throw new IllegalArgumentException( "spearmansRankCorrelationTest: "
  + "sampleSize <= 3" );
//!!!
List< double[2] > ranks = samples;
for ( int i = 0; i < 2; ++i )
{
Collections.sort( ranks, CompareArrays2( i ) );
int j = 1;
while ( j < sampleSize )
{
if ( ranks[ j - 1 ][ i ] != ranks[ j ][ i ] )
{
ranks[ j - 1 ][ i ] = j;
++j;
}
else
{
int k = j + 1;
while ( (k <= sampleSize)
&& (ranks[ k - 1 ][ i ] == ranks[ j - 1 ][ i ]) )
++k;
double aveRank = 0.5 * (j + k - 1);
for ( int m = j; m <= k - 1; ++m )
ranks[ m - 1 ][ i ] = aveRank;
j = k;
}
}
if ( j == sampleSize )
ranks[ j - 1 ][ i ] = j;
}
return LinearCorrelationTest( ranks, pSpearmansR, pT );
*/
        return 0.; //!!!
    }

//-----------------------------------------------------------------------------

    public static
    NormalTestResult
    kendallsTauTest( final List< double[] > samples )
    {
        int sampleSize = samples.size();
        if ( sampleSize <= 1 )
            throw new IllegalArgumentException( "kendallsTauTest: "
                                                + "sampleSize <= 1" );
        int concordantMinusDiscordant = 0;
        int total0 = 0;
        int total1 = 0;
        for ( double[] da1 : samples )
        {
            if ( da1.length != 2 )
                throw new IllegalArgumentException( "kendallsTauTest: "
                                                    + "samples[].size != 2" );
            for ( double[] da2 : samples )
            {
                double diff0 = da1[0] - da2[0];
                double diff1 = da1[1] - da2[1];
                double rel01 = diff0 * diff1;
                if ( rel01 > 0. )   //concordant
                {
                    ++concordantMinusDiscordant;
                    ++total0;
                    ++total1;
                }
                else if ( rel01 < 0. )  //discordant
                {
                    --concordantMinusDiscordant;
                    ++total0;
                    ++total1;
                }
                else //a tie
                {
                    if ( diff0 != 0. )
                        ++total0;
                    if ( diff1 != 0. )
                        ++total1;
                }
            }
        }
        double total01 = total0 * total1;
        assert( total01 > 0. );
        double kendallsTau = concordantMinusDiscordant / Math.sqrt( total01 );
        double variance = (4. * sampleSize  +  10.)
                / (9. * sampleSize * (sampleSize - 1.));
        double prob = ProbDist.normal_DF( kendallsTau, 0., variance );
        //Two-tailed test
        if ( kendallsTau < 0 )
            prob = 2. * prob;
        else
            prob = 2. * (1. - prob);
        return new NormalTestResult( prob, kendallsTau, 0., variance );
    }

//=============================================================================

    public static
    ContingencyTableResult
    chiSquareContingencyTableTest( final Integer[][] table,
                                   boolean yatesCorrection )
    {
        int numRows = table.length;
        if ( numRows <= 1 )
            throw new IllegalArgumentException( "chiSquareContingencyTableTest:"
                                                + " numRows <= 1" );
        int numColumns = table[ 0 ].length;
        if ( numColumns <= 1 )
            throw new IllegalArgumentException( "chiSquareContingencyTableTest:"
                                                + " numColumns <= 1" );
        int sampleTotal = 0;
        int[] rowTotals = new int[ numRows ];
        int[] columnTotals = new int[ numColumns ];
        for ( int i = 0; i < numRows; ++i )
        {
            if ( table[ i ].length != numColumns )
                throw new IllegalArgumentException(
                    "chiSquareContingencyTableTest: row " + i
                    + " size != " + numColumns );
            for ( int j = 0; j < numColumns; ++j )
            {
                int n = table[ i ][ j ];
                sampleTotal += n;
                rowTotals[ i ] += n;
                columnTotals[ j ] += n;
            }
        }
        if ( sampleTotal <= 0 )
            throw new IllegalArgumentException(
                "chiSquareContingencyTableTest: sampleTotal <= 0" );
        int nRows = numRows;
        int nCols = numColumns;
        double chiSquare = 0.;
        double minExpectedCellFreq = sampleTotal;
        for ( int i = 0; i < numRows; ++i )
        {
            if ( rowTotals[ i ] == 0 )
            {
                --nRows;
                continue;
            }
            for ( int j = 0; j < numColumns; ++j )
            {
                if ( columnTotals[ j ] == 0 )
                    continue;
                double expected = (double)rowTotals[i] * columnTotals[j]
                        / sampleTotal;
                if ( expected < minExpectedCellFreq )
                    minExpectedCellFreq = expected;
                if ( expected > 0. )
                {
                    double diff = table[i][j] - expected;
                    if ( yatesCorrection )
                        diff = Math.abs( diff ) - 0.5;
                    chiSquare += diff * diff / expected;
                }
            }
        }
        for ( int j = 0; j < numColumns; ++j )
            if ( columnTotals[ j ] == 0 )
                --nCols;
        assert( nRows > 1 );
        assert( nCols > 1 );
        int degreesOfFreedom = (nRows - 1) * (nCols - 1);

        if ( degreesOfFreedom <= 0 )
            return new ContingencyTableResult(
                new ChiSquareTestResult( 1., chiSquare, degreesOfFreedom ),
                minExpectedCellFreq, sampleTotal );
                
        double prob = 1. - ProbDist.chiSquare_DF( chiSquare, degreesOfFreedom );
        return new ContingencyTableResult(
            new ChiSquareTestResult( prob, chiSquare, degreesOfFreedom ),
            minExpectedCellFreq, sampleTotal );
    }

//.............................................................................

    public static
    ContingencyTableResult
    chiSquareContingencyTableTest( final List< List< Integer > > table,
                                   boolean yatesCorrection )
    {
        int numRows = table.size();
        Integer[][] tableArray = new Integer[ numRows ][];
        for ( int i = 0; i < numRows; ++i )
        {
            List< Integer > row = table.get( i );
            tableArray[ i ] = new Integer[ row.size() ];
            row.toArray( tableArray[ i ] );
        }
        return chiSquareContingencyTableTest( tableArray, yatesCorrection );
    }

//-----------------------------------------------------------------------------

    public static
    double
    fishersExactTest( final Integer[][] table, int tails )
    {
        if ( (table.length != 2) || (table[0].length != 2)
             || (table[1].length != 2) )
            throw new IllegalArgumentException( "fishersExactTest: "
                                                + "table not 2x2" );
        //locate smallest cell value
        int i0 = 0;
        int j0 = 0;
        int minVal = table[0][0];
        for ( int i = 0; i < 2; ++i )
            for ( int j = 0; j < 2; ++j )
                if ( table[i][j] < minVal )
                {
                    i0 = i;
                    j0 = j;
                    minVal = table[i][j];
                }
        int i1 = 1 - i0;
        int j1 = 1 - j0;
        int a = table[i0][j0];
        int b = table[i0][j1];
        int c = table[i1][j0];
        int d = table[i1][j1];
        double prob = ProbDist.hypergeometric_DF( a, (a + b + c + d),
                                                  (a + b), (a + c) );
        if ( tails == 2 )
        {
            if ( ((a + b) == (c + d)) || ((a + c) == (b + d)) )
            {   //symmetrical
                prob *= 2.;
            }
            else
            {
                double prA = ProbDist.hypergeometric_PDF( a, (a + b + c + d),
                                                          (a + b), (a + c) );
                int maxA = Math.min( (a + b), (a + c) );
                for ( int i = maxA; i > a; --i )
                {
                    double pr = ProbDist.hypergeometric_PDF( i, (a + b + c + d),
                                                             (a + b), (a + c) );
                    if ( pr <= prA )
                        prob += pr;
                    else
                        break;
                }
            }
        }
        return prob;
    }

//=============================================================================

    public static
    class NormalTestResult
    {                                                        //NormalTestResult
    //-------------------------------------------------------------------------

        public
        NormalTestResult( double probability, double statistic,
                          double mean, double variance )
        {
            this.probability = probability;
            this.statistic = statistic;
            this.mean = mean;
            this.variance = variance;
        }

    //=========================================================================

        public final double probability;
        public final double statistic;
        public final double mean;
        public final double variance;
        
    //-------------------------------------------------------------------------
    }                                                        //NormalTestResult

//=============================================================================

    public static
    class TTestResult
    {                                                             //TTestResult
    //-------------------------------------------------------------------------

        public
        TTestResult( double probability, double t, int degreesOfFreedom )
        {
            this.probability = probability;
            this.t = t;
            this.degreesOfFreedom = degreesOfFreedom;
        }

    //=========================================================================

        public final double probability;
        public final double t;
        public final int degreesOfFreedom;
        
    //-------------------------------------------------------------------------
    }                                                             //TTestResult

//=============================================================================

    public static
    class ChiSquareTestResult
    {                                                     //ChiSquareTestResult
    //-------------------------------------------------------------------------

        public
        ChiSquareTestResult( double probability,
                             double chiSquare, int degreesOfFreedom )
        {
            this.probability = probability;
            this.chiSquare = chiSquare;
            this.degreesOfFreedom = degreesOfFreedom;
        }

    //=========================================================================

        public final double probability;
        public final double chiSquare;
        public final int degreesOfFreedom;
        
    //-------------------------------------------------------------------------
    }                                                     //ChiSquareTestResult

//=============================================================================

    public static
    class FTestResult
    {                                                             //FTestResult
    //-------------------------------------------------------------------------

        public
        FTestResult( double probability, double f,
                     int degreesOfFreedom1, int degreesOfFreedom2 )
        {
            this.probability = probability;
            this.f = f;
            this.degreesOfFreedom1 = degreesOfFreedom1;
            this.degreesOfFreedom2 = degreesOfFreedom2;
        }

    //=========================================================================

        public final double probability;
        public final double f;
        public final int degreesOfFreedom1;
        public final int degreesOfFreedom2;
        
    //-------------------------------------------------------------------------
    }                                                             //FTestResult

//=============================================================================

    public static
    class CorrelationTestResult
    {                                                   //CorrelationTestResult
    //-------------------------------------------------------------------------

        public
        CorrelationTestResult( TTestResult tResult,
                               double r,
                               double mean0, double mean1,
                               double variance0, double variance1 )
        {
            this.tResult = tResult;
            this.r = r;
            this.mean0 = mean0;
            this.mean1 = mean1;
            this.variance0 = variance0;
            this.variance1 = variance1;
        }

    //=========================================================================

        public final TTestResult tResult;
        public final double r;
        public final double mean0;
        public final double mean1;
        public final double variance0;
        public final double variance1;
        
    //-------------------------------------------------------------------------
    }                                                   //CorrelationTestResult

//=============================================================================

    public static
    class SimpleRegressionResult
    {                                                  //SimpleRegressionResult
    //-------------------------------------------------------------------------

        public
        SimpleRegressionResult( double alpha, double beta,
                                double alphaVariance, double betaVariance,
                                double residualVariance )
        {
            this.alpha = alpha;
            this.beta = beta;
            this.alphaVariance = alphaVariance;
            this.betaVariance = betaVariance;
            this.residualVariance = residualVariance;
        }
        
    //=========================================================================

        public final double alpha;
        public final double beta;
        public final double alphaVariance;
        public final double betaVariance;
        public final double residualVariance;
        
    //-------------------------------------------------------------------------
    }                                                  //SimpleRegressionResult

//=============================================================================

    public static
    class ContingencyTableResult
    {                                                  //ContingencyTableResult
    //-------------------------------------------------------------------------

        public
        ContingencyTableResult( ChiSquareTestResult chiSquareResult,
                                double minExpectedCellFreq,
                                int sampleTotal )
        {
            this.chiSquareResult = chiSquareResult;
            this.minExpectedCellFreq = minExpectedCellFreq;
            this.sampleTotal = sampleTotal;
        }

    //=========================================================================

        public final ChiSquareTestResult chiSquareResult;
        public final double minExpectedCellFreq;
        public final int sampleTotal;
        
    //-------------------------------------------------------------------------
    }                                                  //ContingencyTableResult


//-----------------------------------------------------------------------------
}                                                                       //Stats


//*****************************************************************************
