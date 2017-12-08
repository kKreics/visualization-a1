/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import static java.lang.Math.abs;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private int interactiveStep = 35;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset
        tFunc.setTestFunc();


        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }

    VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] >= gradients.getDimX()
                || coord[1] < 0 || coord[1] >= gradients.getDimY()
                || coord[2] < 0 || coord[2] >= gradients.getDimZ()) {
            return new VoxelGradient();
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return gradients.getGradient(x, y, z);
    }


    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] + 1 >= volume.getDimX() || coord[1] < 0 || coord[1] + 1 >= volume.getDimY()
                || coord[2] < 0 || coord[2] + 1 >= volume.getDimZ()) {
            return -1;
        }

        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);

        if (panel.shadingActive) return volume.getVoxel(x0, y0, z0);

        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        double alfa = coord[0] - x0;
        double beta = coord[1] - y0;
        double gamma = coord[2] - z0;


        return (short) (
            ((1 - alfa) * (1 - beta) * (1 - gamma) * volume.getVoxel(x0, y0, z0)) +
            (alfa * (1 - beta) * (1 - gamma) * volume.getVoxel(x1, y0, z0)) +
            ((1 - alfa) * beta * (1 - gamma) * volume.getVoxel(x0, y1, z0)) +
            (alfa * beta * (1 - gamma) * volume.getVoxel(x1, y1, z0)) +
            ((1 - alfa) * (1 - beta) * gamma * volume.getVoxel(x0, y0, z1)) +
            (alfa * (1 - beta) * gamma * volume.getVoxel(x1, y0, z1)) +
            ((1 - alfa) * beta * gamma * volume.getVoxel(x0, y1, z1)) +
            (alfa * beta * gamma * volume.getVoxel(x1, y1, z1))
        );
    }


    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();


        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    /*
    Calculate shortest distance from point "P" to the line determined
    by points "A" and "B".

    See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    */
    double GetDistance(double[] P, double[] A, double[] B) {
        double[] PMinusA = new double[3];
        double[] PMinusB = new double[3];
        double[] BMinusA = new double[3];

        VectorMath.setVector(PMinusA, P[0] - A[0], P[1] - A[1], P[2] - A[2]);
        VectorMath.setVector(PMinusB, P[0] - B[0], P[1] - B[1], P[2] - B[2]);
        VectorMath.setVector(BMinusA, B[0] - A[0], B[1] - A[1], B[2] - A[2]);

        double[] crossProd = new double[3];
        crossProd = VectorMath.crossproduct(PMinusA, PMinusB, crossProd);

        double numerator = VectorMath.length(crossProd);
        double denominator = VectorMath.length(BMinusA);

        return numerator/denominator;
    }


    void mip(double[] viewMatrix) {
        // clear image
        int step = interactiveMode ? interactiveStep : 1;
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();


        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                short currentMaxIntensity = 0;

                // Move in the direction of viewVec, starting in the plane
                // that cuts the figure in two, perpendicular to the viewPlane
                double[] pixelCoordViewVec = VectorMath.cloneVector(pixelCoord);
                double[] viewScale = VectorMath.scaleVector(viewVec, step);

                while (true) {
                    short voxelIntensity = getVoxel(pixelCoordViewVec);

                    // We are out of the figure
                    if (voxelIntensity == -1) {
                        break;
                    }

                    if (voxelIntensity > currentMaxIntensity) {
                        currentMaxIntensity = voxelIntensity;
                    }

                    pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);
                }

                if (panel.enableRaycastBothWays) {
                    // Move in the opposite direction of viewVec, starting in the
                    // plane that cuts the figure in two, perpendicular to the
                    // viewPlane
                    pixelCoordViewVec = VectorMath.cloneVector(pixelCoord);
                    viewScale = VectorMath.scaleVector(viewVec, -step);

                    while (true) {
                        // When moving in the opposite direction, first apply the
                        // increment so as not to do calculations with the middle
                        // pixel twice
                        pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);

                        short voxelIntensity = getVoxel(pixelCoordViewVec);

                        // We are out of the figure
                        if (voxelIntensity == -1) {
                            break;
                        }

                        if (voxelIntensity > currentMaxIntensity) {
                            currentMaxIntensity = voxelIntensity;
                        }
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = currentMaxIntensity/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = currentMaxIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
//                voxelColor = tFunc.getColor(currentMaxIntensity);


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void composite(double[] viewMatrix) {
        // clear image
        int step = interactiveMode ? interactiveStep : 1;
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                TFColor voxelColor = new TFColor();
                TFColor tmpColor = new TFColor();

                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                // Move in the direction of viewVec, starting in the plane
                // that cuts the figure in two, perpendicular to the viewPlane
                double[] pixelCoordViewVec = VectorMath.cloneVector(pixelCoord);
                double[] viewScale = VectorMath.scaleVector(viewVec, step);

                // Advance until the end of the volume (out of bounds)
                while (getVoxel(pixelCoordViewVec) != -1) {
                    pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);
                }

                // Reverse the moving vector
                viewScale = VectorMath.scaleVector(viewScale, -1);
                // Get inside the volume again
                pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);

                // Iterate the volume until we get out of the volume from the other side
                while (true) {
                    short voxelIntensity = getVoxel(pixelCoordViewVec);

                    // We are out of the figure
                    if (voxelIntensity == -1) {
                        break;
                    }

                    tmpColor = tFunc.getColor(voxelIntensity);
                    voxelColor.r = tmpColor.r * tmpColor.a + (1 - tmpColor.a) * voxelColor.r;
                    voxelColor.g = tmpColor.g * tmpColor.a + (1 - tmpColor.a) * voxelColor.g;
                    voxelColor.b = tmpColor.b * tmpColor.a + (1 - tmpColor.a) * voxelColor.b;

                    pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);
                }

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    // Implementation of the Levoy's formula from the paper. Returns the
    // opacity value for some voxel.
    //
    // Parameters
    // ----------
    // coord
    //     Coordinates of the point for which to apply the formula. I think in
    //     the paper they are represented by x_i.
    double levoysFormula(double[] coord) {
        // Specifies the value for the voxels which we want to focus on
        int f_v = tfEditor2D.triangleWidget.baseIntensity;
        // Desired thickness in the voxels
        double r = tfEditor2D.triangleWidget.radius;
        // Voxel value of the specified point
        short f_xi = getVoxel(coord);

        // Gradient value of the specified point
        // Note that its magnitude (mag) is always positive, so we do not have
        // to put so many abs(...) around it like in the book
        VoxelGradient gradient = getGradient(coord);

        boolean condition_one = (gradient.mag == 0) && (f_xi == f_v);
        boolean condition_two = (gradient.mag > 0)
                                && (f_xi - r*gradient.mag <= f_v)
                                && (f_xi + r*gradient.mag >= f_v);
        if (gradient.mag < tfEditor2D.triangleWidget.minGradientValue
                || gradient.mag > tfEditor2D.triangleWidget.maxGradientValue) {
            return 0;
        } else if (condition_one) {
            return 1;
        }
        else if (condition_two) {
            double division_part = (f_v - f_xi) / gradient.mag;

            return 1 - (1/r) * abs(division_part);
        }
        else {
            return 0;
        }
    }

    // l and v are viewVecs, g is gradient
    TFColor phong(double[] v, double[] g, TFColor color) {
        // take the parameters from the assignment description
        double kAmb = 0.1;
        double kDiff = 0.7;
        double kSpec = 0.2;
        double alpha = 10;

        // dot product here always returns a non negative value
        double x = Math.max(0, VectorMath.dotproduct(v, g));
        double y = Math.pow(x, alpha);

        // use simplified model formula from slide
        return new TFColor(
            kAmb * color.r + kDiff * color.r * x + kSpec * y,
            kAmb * color.g + kDiff * color.g * x + kSpec * y,
            kAmb * color.b + kDiff * color.b * x + kSpec * y,
            color.a
        );
    }


    void gradientOpacity(double[] viewMatrix) {
        // clear image
        int step = interactiveMode ? 12 : 1;
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        viewVec = VectorMath.normalized(viewVec);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                TFColor voxelColor = new TFColor();
                TFColor tmpColor = new TFColor();
                //TFColor shadeColor = new TFColor();

                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                // Move in the direction of viewVec, starting in the plane
                // that cuts the figure in two, perpendicular to the viewPlane
                double[] pixelCoordViewVec = VectorMath.cloneVector(pixelCoord);
                double[] viewScale = VectorMath.scaleVector(viewVec, step);
                viewScale = VectorMath.scaleVector(viewScale, -1);

                // Advance until the end of the volume (out of bounds)
                while (getVoxel(pixelCoordViewVec) != -1) {
                    pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);
                }

                // Reverse the moving vector
                viewScale = VectorMath.scaleVector(viewScale, -1);
                // Get inside the volume again
                pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);

                // Iterate the volume until we get out of the volume from the other side
                voxelColor.r = 0;
                voxelColor.b = 0;
                voxelColor.g = 0;
                voxelColor.a = 1;
                while (true) {
                    short voxelIntensity = getVoxel(pixelCoordViewVec);

                    // We are out of the figure
                    if (voxelIntensity == -1) {
                        break;
                    }

                    tmpColor = tfEditor2D.triangleWidget.color;
                    tmpColor.a = levoysFormula(pixelCoordViewVec);

                    if (panel.shadingActive) {
                        VoxelGradient h = getGradient(pixelCoordViewVec);

                        if (h.mag != 0) {
                            tmpColor = phong(viewVec, h.normal(), tmpColor);
                        }
                    }

                    voxelColor.r = tmpColor.r * tmpColor.a + (1 - tmpColor.a) * voxelColor.r;
                    voxelColor.g = tmpColor.g * tmpColor.a + (1 - tmpColor.a) * voxelColor.g;
                    voxelColor.b = tmpColor.b * tmpColor.a + (1 - tmpColor.a) * voxelColor.b;

                    pixelCoordViewVec = VectorMath.addVectors(pixelCoordViewVec, viewScale);
                }

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        // Here we should change the "mip" call to whichever visualization
        // method we want (i.e.: mip, slicer...).
        switch (panel.selectedOption) {
            case 0:
                slicer(viewMatrix);
                break;
            case 1:
                mip(viewMatrix);
                break;
            case 2:
                composite(viewMatrix);
                break;
            case 3:
                gradientOpacity(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
        }

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
