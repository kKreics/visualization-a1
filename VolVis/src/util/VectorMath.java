/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package util;

/**
 *
 * @author michel
 */
public class VectorMath {

    // assign coefficients c0..c2 to vector v
    public static void setVector(double[] v, double c0, double c1, double c2) {
        v[0] = c0;
        v[1] = c1;
        v[2] = c2;
    }

    // Return the sum of two vectors, element-wise
    public static double[] addVectors(double[] v, double[] w) {
        double[] result = new double[3];
        VectorMath.setVector(result, v[0] + w[0], v[1] + w[1], v[2] + w[2]);

        return result;
    }

    // compute dotproduct of vectors v and w
    public static double dotproduct(double[] v, double[] w) {
        double r = 0;
        for (int i=0; i<3; i++) {
            r += v[i] * w[i];
        }
        return r;
    }

     // Return a clone of the given vector
    public static double[] cloneVector(double[] v) {
        double[] result = new double[3];
        VectorMath.setVector(result, v[0], v[1], v[2]);

        return result;
    }

    // Return a scaled version of a given vector by some factor
    public static double[] scaleVector(double[] v, double factor) {
        double[] result = new double[3];
        VectorMath.setVector(result, v[0]*factor, v[1]*factor, v[2]*factor);

        return result;
    }

    // compute distance between vectors v and w
    public static double distance(double[] v, double[] w) {
        double[] tmp = new double[3];
        VectorMath.setVector(tmp, v[0]-w[0], v[1]-w[1], v[2]-w[2]);
        return Math.sqrt(VectorMath.dotproduct(tmp, tmp));
    }

    // compute dotproduct of v and w
    public static double[] crossproduct(double[] v, double[] w, double[] r) {
        r[0] = v[1] * w[2] - v[2] * w[1];
        r[1] = v[2] * w[0] - v[0] * w[2];
        r[2] = v[0] * w[1] - v[1] * w[0];
        return r;
    }
    
    // compute length of vector v
    public static double length(double[] v) {
        return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    

}
