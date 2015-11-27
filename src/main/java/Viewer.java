/* This code is adapted and mofified from PTViewer.
 * The original PTViewer code is copyrighted by Helmul Dersch.
 * PTViewer Copyright (C) 2000  - Helmut Dersch  der@fh-furtwangen.de
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  
 */

import ij.IJ;
import ij.plugin.frame.PlugInFrame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.MemoryImageSource;

public class Viewer extends PlugInFrame implements MouseListener, MouseMotionListener, KeyListener {

    private static final String title = "Equirectangular Viewer";

    /**
     * Minimum horizontal field-of-view
     */
    public static final double HFOV_MIN = 10.5;

    /**
     * Maximum horizontal field-of-view
     */
    public static final double HFOV_MAX = 165.0;

    /**
     * Panorama Viewport
     */
    private Image view = null;

    /**
     * Offscreen image
     */
    private Image offImage = null;

    /**
     * Offscreen image's rendering graphics
     */
    private Graphics offGraphics = null;

    /**
     * Offscreen image dimension
     */
    private int offwidth = 0;
    private int offheight = 0;

    /**
     * Memory image source
     */
    MemoryImageSource source = null;

    /**
     * Viewer window dimension
     */
    private int vwidth = 0;
    private int vheight = 0;

    /**
     * Position of viewer window
     */
    private int vx = 0;
    private int vy = 0;

    /**
     * Dimension of equirectangular panorama image
     */
    private int pwidth = 0;
    private int pheight = 0;

    /**
     * RGB data of viewer image
     */
    private int[] vdata = null;

    /**
     * RGB data of equirectangular panorama image
     */
    private int[][] pdata = null;

    /**
     * Initialize flag
     */
    private boolean inited = false;

    /**
     * Current yaw angle
     */
    public double yaw = 0.0;

    /**
     * Current horizontal field of view
     */
    public double hfov = 70.0;

    /**
     * Maximum horizontal field of view
     */
    public double hfov_min = HFOV_MIN;

    /**
     * Minimum horizontal field of view
     */
    public double hfov_max = HFOV_MAX;

    /**
     * Current tilt angle
     */
    public double pitch = 0.0;

    /**
     * Maximum tilt angle
     */
    public double pitch_max = 90.0;

    /**
     * Minimum tilt angle
     */
    public double pitch_min = -90.0;

    /**
     * Maximum pan angle
     */
    public double yaw_max = 180.0;

    /**
     * Minimum pan angle
     */
    public double yaw_min = -180.0;

    /**
     * Current zoom factor
     */
    public double zoom = 1.0;

    /**
     * Number of steps for panning one hfov
     */
    public double pan_steps = 20.0;

    /**
     * Mousebutton has been pressed
     */
    private boolean panning = false;

    /**
     * Image needs recalculation
     */
    public boolean dirty = true;

    /**
     * Mouse coordinates
     */
    private int oldx = 0, oldy = 0;
    private int newx = 0, newy = 0;

    /**
     * Size of Lookup table for atan2 routines = 2^12
     */
    private static final int NATAN = 4096;

    /**
     * Size of Lookup table for sqrt routine = 2^12
     */
    private static final int NSQRT = 4096;

    /**
     * 3 x 3 Transformation Matrix
     */
    private double mt[][] = null;

    /**
     * Integer version of matrix above
     */
    private int mi[][] = null;

    /**
     * Atan function lookup table
     */
    private int atan_LU_HR[] = null;

    /**
     * Square root lookup table
     */
    private int sqrt_LU[] = null;

    /**
     * Bilinear interpolation weights
     */
    private int mweights[][] = null;

    /**
     * atan2 special value with x = 0
     */
    private int PV_atan0_HR;

    /**
     * atan2 special value with y = 0
     */
    private int PV_pi_HR;

    final static private String NL = System.getProperty("line.separator");

    public Viewer(int[][] pd) {
        super(title);
        pdata = pd;
        pheight = pdata.length;
        pwidth = pdata[0].length;
    }

    public void init() {
        math_setLookUp(pdata);
        math_init();
        addMouseListener(this);
        addMouseMotionListener(this);
        addKeyListener(this);
        inited = true;
    }

    @Override
    public void update(Graphics g) {
        paint(g);
    }

    @Override
    public void paint(Graphics g) {
        if (!inited) {
            return;
        }

        if (offImage == null) {
            if (offwidth == 0 || offheight == 0) {
                offwidth = getSize().width;
                offheight = getSize().height;
            }
            offImage = createImage(offwidth, offheight);
            offGraphics = offImage.getGraphics();
        }

        if (vdata == null) {
            if (vwidth == 0 || vheight == 0) {
                vwidth = getSize().width;
                vheight = getSize().height;
            }

            while (math_fovy(hfov, vwidth, vheight) > pitch_max - pitch_min) {
                hfov /= 1.03;
            }

            double fovy2 = math_fovy(hfov, vwidth, vheight) / 2.0;

            if (pitch > pitch_max - fovy2 && pitch_max != 90.0) {
                pitch = 0.0;
            }

            if (pitch < pitch_min + fovy2 && pitch_min != -90.0) {
                pitch = 0.0;
            }

            vdata = new int[vwidth * vheight];

            dirty = true;

            source = new MemoryImageSource(vwidth, vheight, vdata, 0, vwidth);
            source.setAnimated(true);

            if (view == null) {
                view = createImage(source);
            }
        }

        if (panning) {
            double scale = 1.0 / 2000.0 * hfov / 70.0 * 320.0 / vwidth;
            gotoView(yaw + scale * (double) ((newx - oldx) * (newx - oldx)) * (newx > oldx ? 1.0 : -1.0),
                    pitch + scale * (double) ((oldy - newy) * (oldy - newy)) * (oldy > newy ? 1.0 : -1.0),
                    hfov * zoom);
        }

        if (dirty) {
            for (int i = 0; i < vdata.length; i++) {
                vdata[i] = 0;
            }

            math_extractview(pdata, vdata, vwidth, hfov, yaw, pitch);

            dirty = false;
            source.newPixels();
        }

        offGraphics.drawImage(view, vx, vy, this);
        g.drawImage(offImage, 0, 0, this);
    }

    /**
     * Go to a new view
     *
     * @param pan pan angle
     * @param tilt tilt angle
     * @param fov field-of-view
     */
    public void gotoView(double pan, double tilt, double fov) {
        if (pan == yaw && tilt == pitch && fov == hfov) {
            return;
        }

        while (pan > 180.0) {
            pan -= 360.0;
        }

        while (pan < -180.0) {
            pan += 360.0;
        }

        final double f = math_fovy(fov, vwidth, vheight) / 2.0;

        if (tilt > pitch_max - f && pitch_max != 90.0) {
            tilt = pitch_max - f;
        } else if (tilt > pitch_max) {
            tilt = pitch_max;
        } else if (tilt < pitch_min + f && pitch_min != -90.0) {
            tilt = pitch_min + f;
        } else if (tilt < pitch_min) {
            tilt = pitch_min;
        }

        // Check and correct for yaw_max/min
        if (yaw_max != 180.0 || yaw_min != -180.0) {
            double x[];
            double xl, xr;

            // check left edge
            x = math_view2pano(0, (pitch > 0.0 ? 0 : vheight - 1), vwidth, vheight,
                    pwidth, pheight, pan, tilt, fov);
            xl = x[0];

            x = math_view2pano(vwidth - 1, (pitch > 0.0 ? 0 : vheight - 1), vwidth, vheight,
                    pwidth, pheight, pan, tilt, fov);
            xr = x[0];

            if ((xr - xl) > (yaw_max - yaw_min) / 360.0 * (double) pwidth) {
                return;
            }

            if (xl < (yaw_min + 180.0) / 360.0 * (double) pwidth) {
                pan += yaw_min - (xl / (double) pwidth * 360.0 - 180.0);
            }

            if (xr > (yaw_max + 180.0) / 360.0 * (double) pwidth) {

                pan -= (xr / (double) pwidth * 360.0 - 180.0) - yaw_max;
            }
        }

        if (2.0 * f <= pitch_max - pitch_min
                && fov <= hfov_max
                && fov >= hfov_min
                && fov <= yaw_max - yaw_min
                && tilt <= pitch_max
                && tilt >= pitch_min
                && pan <= yaw_max
                && pan >= yaw_min) {

            if (pan != yaw || tilt != pitch || fov != hfov) {
                yaw = pan;
                pitch = tilt;
                hfov = fov;
                dirty = true;
                repaint();
            }
        }
    }

    private void math_init() {
        mt = new double[3][3];
        mi = new int[3][3];
        // Set up Multiplication table
        mweights = new int[256][256];
        for (int i = 0; i < 256; i++) {
            for (int k = 0; k < 256; k++) {
                mweights[i][k] = i * k;
            }
        }
    }

    private void math_setLookUp(int[][] pd) {
        double z, dz;

        if (pd == null) {
            return;
        }

        final int ph = pd.length;
        final int pw = pd[0].length;

        if (atan_LU_HR == null) {
            atan_LU_HR = new int[NATAN + 1];
            sqrt_LU = new int[NSQRT + 1];

            dz = 1.0 / (double) NSQRT;
            z = 0.0;

            for (int i = 0; i < NSQRT; i++, z += dz) {
                sqrt_LU[i] = (int) (Math.sqrt(1.0 + z * z) * NSQRT);
            }

            sqrt_LU[NSQRT] = (int) (Math.sqrt(2.0) * NSQRT);

        }

        dz = 1.0 / (double) NATAN;
        z = 0.0;

        final double dist_e = (double) pw / (2.0 * Math.PI);

        PV_atan0_HR = pw * 64;
        PV_pi_HR = pw * 128;

        for (int i = 0; i < NATAN + 1; i++, z += dz) {
            if (i < NATAN) {
                atan_LU_HR[i] = (int) (dist_e * Math.atan(z / (1.0 - z)) * 256.0 + 0.5);
            } else {
                atan_LU_HR[i] = (int) (dist_e * Math.PI / 2.0 * 256.0 + 0.5);
            }
        }
    }

    /**
     * Compute vertical field-of-view from horizontal field-of-view
     *
     * @param Hfov horizontal field-of-view
     * @param vw viewer width
     * @param vh viewer height
     * @return vertical field-of-view
     */
    private double math_fovy(double Hfov, int vw, int vh) {
        return (360.0 / Math.PI) * Math.atan((double) vh / (double) vw * Math.tan(Hfov / 2.0 * Math.PI / 180.0));
    }

    /**
     * Square root computation using lookup table
     *
     * @param x1
     * @param x2
     * @return square root
     */
    private int my_sqrt(int x1, int x2) {
        return x1 > x2 ? (x1 * sqrt_LU[(x2 << 12) / x1]) >> 12 : x2 == 0 ? 0 : (x2 * sqrt_LU[(x1 << 12) / x2]) >> 12;
    }

    /**
     * Atan2 computation usin lookup table
     *
     * @param y
     * @param x
     * @return atan2 value
     */
    private int my_atan2_HR(int y, int x) {
        if (x > 0) {
            if (y > 0) {
                return atan_LU_HR[(int) (NATAN * y / (x + y))];
            } else {
                return -atan_LU_HR[(int) (NATAN * (-y) / (x - y))];
            }
        }

        if (x == 0) {
            if (y > 0) {
                return PV_atan0_HR;
            } else {
                return -PV_atan0_HR;
            }
        }

        if (y < 0) {
            return atan_LU_HR[(int) (NATAN * y / (x + y))] - PV_pi_HR;
        } else {
            return -atan_LU_HR[(int) (NATAN * (-y) / (x - y))] + PV_pi_HR;
        }
    }

    /**
     * Bilinear interpolation
     *
     * @param p00 pixel p00
     * @param p01 pixel p01
     * @param p10 pixel p10
     * @param p11 pixel p11
     * @param dx fraction in x direction
     * @param dy fraction in y direction
     * @return interpolated pixel value
     */
    final int bilinearInterpolation(int p00, int p01, int p10, int p11, int dx, int dy) {

        int w1[], w2[], rd, gd, bd, yr, yg, yb, weight;

        w1 = mweights[dx & 0xff];
        w2 = mweights[(255 - dx) & 0xff];

        rd = w2[(p00 >> 16) & 0xff] + w1[(p01 >> 16) & 0xff];
        gd = w2[(p00 >> 8) & 0xff] + w1[(p01 >> 8) & 0xff];
        bd = w2[(p00) & 0xff] + w1[(p01) & 0xff];

        yr = w2[(p10 >> 16) & 0xff] + w1[(p11 >> 16) & 0xff];
        yg = w2[(p10 >> 8) & 0xff] + w1[(p11 >> 8) & 0xff];
        yb = w2[(p10) & 0xff] + w1[(p11) & 0xff];

        weight = 255 - dy;

        rd = (rd * weight + yr * dy) & 0xff0000;
        gd = (gd * weight + yg * dy) >> 16;
        bd = (bd * weight + yb * dy) >> 16;

        return (rd + (gd << 8) + bd + 0xff000000);
    }

    /**
     * Set rotation matrix
     *
     * @param a tilt angle
     * @param b pan angle
     * @param m rotaton matrix
     */
    private void SetMatrix(double a, double b, double m[][]) {
        double mx[][], my[][];

        mx = new double[3][3];
        my = new double[3][3];

        mx[0][0] = 1.0;
        mx[0][1] = 0.0;
        mx[0][2] = 0.0;
        mx[1][0] = 0.0;
        mx[1][1] = Math.cos(a);
        mx[1][2] = Math.sin(a);
        mx[2][0] = 0.0;
        mx[2][1] = -mx[1][2];
        mx[2][2] = mx[1][1];

        my[0][0] = Math.cos(b);
        my[0][1] = 0.0;
        my[0][2] = -Math.sin(b);
        my[1][0] = 0.0;
        my[1][1] = 1.0;
        my[1][2] = 0.0;
        my[2][0] = -my[0][2];
        my[2][1] = 0.0;
        my[2][2] = my[0][0];

        matrix_matrix_mult(mx, my, m);
    }

    /**
     * Matrix multiply
     *
     * @param m1 rotation matrix
     * @param m2 rotation matrix
     * @param result composition rotation matrix
     */
    final void matrix_matrix_mult(double m1[][], double m2[][], double result[][]) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                result[i][k] = m1[i][0] * m2[0][k] + m1[i][1] * m2[1][k] + m1[i][2] * m2[2][k];
            }
        }
    }

    /**
     * Extract a new view based on pan, tilt and fov angle.
     *
     * @param pd panorama data
     * @param vd viewer data
     * @param vw viewer width
     * @param fov field-of-view
     * @param pan pan angle
     * @param tilt titl angle
     */
    final void math_extractview(int[][] pd, int[] vd, int vw, double fov, double pan, double tilt) {
        math_set_int_matrix(fov, pan, tilt, vw);
        math_transform(pd, pd[0].length, pd.length, vd, vw, vd.length / vw);
    }

    /**
     * Set rotaton maxtrix
     *
     * @param fov field-of-view
     * @param pan pan angle
     * @param tilt titl angle
     * @param vw viewer width
     */
    final void math_set_int_matrix(double fov, double pan, double tilt, int vw) {
        double a = fov * 2.0 * Math.PI / 360.0;
        double p = (double) vw / (2.0 * Math.tan(a / 2.0));

        SetMatrix(tilt * 2.0 * Math.PI / 360.0, pan * 2.0 * Math.PI / 360.0, mt);

        mt[0][0] /= p;
        mt[0][1] /= p;
        mt[0][2] /= p;
        mt[1][0] /= p;
        mt[1][1] /= p;
        mt[1][2] /= p;

        double ta = (a > 0.3 ? 65536.0 * 2.0 / a : 65536.0 * 2.0 / 0.3);

        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                mi[i][k] = (int) (ta * mt[i][k] + 0.5);
            }
        }

    }

    /**
     * Inverse transfrom from viewer to equirectangular panorama image
     *
     * @param pd equirectangular panorama image
     * @param pw equirectangular panorama width
     * @param ph equirectangular panorama height
     * @param vd viewer data
     * @param vw viewer width
     * @param vh viewer height
     */
    final void math_transform(int[][] pd, int pw, int ph, int[] vd, int vw, int vh) {
        final int mix = pw - 1;
        final int miy = ph - 1;

        // Variables used to convert screen coordinates to cartesian coordinates
        final int vw2 = (vw - 1) / 2;
        final int vh2 = vh / 2;
        final int sw2 = pw / 2;
        final int sh2 = ph / 2;

        final int x_min = -vw2;
        final int x_max = vw - vw2;
        final int y_min = -vh2;
        final int y_max = vh - vh2;

        // Variables used to perform bilinear interpolation
        int p00, p01, p10, p11;
        int dx, dy;
        int idp;

        int cy = 0;

        // 3D coordinates
        int v0, v1, v2;

        // Increment matrix elements
        int m0 = mi[1][0] * y_min + mi[2][0];
        int m1 = mi[1][1] * y_min + mi[2][1];
        int m2 = mi[1][2] * y_min + mi[2][2];

        int idx;

        // 2D equirectangular image coordinates
        int xs, ys;

        for (int y = y_min; y < y_max; y++, cy += vw, m0 += mi[1][0], m1 += mi[1][1], m2 += mi[1][2]) {
            idx = cy;

            v0 = mi[0][0] * x_min + m0;
            v1 = m1; // mi[0][1] = 0
            v2 = mi[0][2] * x_min + m2;

            for (int x = x_min; x < x_max; x++, idx++, v0 += mi[0][0], v2 += mi[0][2]) {
                if (vd[idx] != 0) {
                    continue;
                }

                xs = my_atan2_HR(v0, v2);
                ys = my_atan2_HR(v1, my_sqrt(Math.abs(v2), Math.abs(v0)));

                dx = xs & 0xff;
                dy = ys & 0xff;

                xs = (xs >> 8) + sw2;
                ys = (ys >> 8) + sh2;

                if (ys >= 0 && ys < miy && xs >= 0 && xs < mix) {
                    p00 = pd[ys][xs];
                    p01 = pd[ys][xs + 1];
                    p10 = pd[ys + 1][xs];
                    p11 = pd[ys + 1][xs + 1];
                } else {
                    if (ys < 0) {
                        idp = 0;
                    } else if (ys > miy) {
                        idp = miy;
                    } else {
                        idp = ys;
                    }

                    if (xs < 0) {
                        p00 = pd[idp][mix];
                    } else if (xs > mix) {
                        p00 = pd[idp][0];
                    } else {
                        p00 = pd[idp][xs];
                    }

                    xs += 1;
                    if (xs < 0) {
                        p01 = pd[idp][mix];
                    } else if (xs > mix) {
                        p01 = pd[idp][0];
                    } else {
                        p01 = pd[idp][xs];
                    }
                    xs -= 1;

                    ys += 1;
                    if (ys < 0) {
                        idp = 0;
                    } else if (ys > miy) {
                        idp = miy;
                    } else {
                        idp = ys;
                    }

                    if (xs < 0) {
                        p10 = pd[idp][mix];
                    } else if (xs > mix) {
                        p10 = pd[idp][0];
                    } else {
                        p10 = pd[idp][xs];
                    }

                    xs += 1;
                    if (xs < 0) {
                        p11 = pd[idp][mix];
                    } else if (xs > mix) {
                        p11 = pd[idp][0];
                    } else {
                        p11 = pd[idp][xs];
                    }
                }
                vd[idx] = bilinearInterpolation(p00, p01, p10, p11, dx, dy);
            }
        }
    }

    final double[] math_view2pano(int xv, int yv, int vw, int vh,
            int pw, int ph,
            double pan, double tilt, double fov) {
        double a, p, dr;
        double v0, v1, v2;

        double dist_e = (double) pw / (2.0 * Math.PI);
        a = fov * 2.0 * Math.PI / 360.0;	// field of view in rad		
        p = (double) vw / (2.0 * Math.tan(a / 2.0));
        dr = (int) (p + .5);

        SetMatrix(tilt * 2.0 * Math.PI / 360.0, pan * 2.0 * Math.PI / 360.0, mt);

        xv -= vw / 2;
        yv -= vh / 2;

        v0 = mt[0][0] * xv + mt[1][0] * yv + mt[2][0] * dr;
        v1 = mt[0][1] * xv + mt[1][1] * yv + mt[2][1] * dr;
        v2 = mt[0][2] * xv + mt[1][2] * yv + mt[2][2] * dr;

        double[] xp = new double[2];

        xp[0] = dist_e * Math.atan2(v0, v2) + pw / 2.0;
        xp[1] = dist_e * Math.atan2(v1, Math.sqrt(v2 * v2 + v0 * v0)) + ph / 2.0;

        return xp;
    }

    @Override
    public void mouseClicked(MouseEvent e) {

    }

    @Override
    public void mousePressed(MouseEvent e) {
        final int x = e.getX();
        final int y = e.getY();

        if (!panning) {
            oldx = x;
            oldy = y;
            panning = true;
            if (e.isShiftDown()) {
                zoom = 1.0 / 1.03;
            } else if (e.isControlDown()) {
                zoom = 1.03;
            } else {
                zoom = 1.0;
            }
            repaint();
        }

        newx = x;
        newy = y;
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        final int x = e.getX();
        final int y = e.getY();

        newx = x;
        newy = y;
        panning = false;
        zoom = 1.0;
    }

    @Override
    public void mouseEntered(MouseEvent e) {

    }

    @Override
    public void mouseExited(MouseEvent e) {

    }

    @Override
    public void mouseDragged(MouseEvent e) {
        final int x = e.getX();
        final int y = e.getY();
        panning = true;
        newx = x;
        newy = y;

        repaint();
    }

    @Override
    public void mouseMoved(MouseEvent e) {

    }

    @Override
    public void keyTyped(KeyEvent e) {

    }

    @Override
    public void keyPressed(KeyEvent e) {
        final int key = e.getKeyCode();

        switch (key) {
            case KeyEvent.VK_UP:
                panUp();
                break;
            case KeyEvent.VK_DOWN:
                panDown();
                break;
            case KeyEvent.VK_LEFT:
                panLeft();
                break;
            case KeyEvent.VK_RIGHT:
                panRight();
                break;
            case KeyEvent.VK_EQUALS:
                zoomIn();
                break;
            case KeyEvent.VK_MINUS:
                zoomOut();
                break;
            case KeyEvent.VK_F1:
            case KeyEvent.VK_H:
                IJ.showMessage(
                        "Interactive Equirectangular Panorama Viewer",
                        "Mouse control:" + NL + " " + NL
                        + "Pan and tilt the panorama by dragging the image in the viewer window and" + NL
                        + "zoom in and out using the Shift and Control key with left mouse key." + NL + " " + NL
                        + "Key control:" + NL + " " + NL
                        + "CURSOR LEFT - Pan left." + NL
                        + "CURSOR RIGHT - Pan right." + NL
                        + "CURSOR UP - Tilt up." + NL
                        + "CURSOR DOWN - Tilt down." + NL
                        + "= key - Zoom in." + NL
                        + "- key - Zoom out.");
                break;
        }
    }

    @Override
    public void keyReleased(KeyEvent e) {

    }

    /**
     * Zoom in 3%
     */
    public void zoomIn() {
        gotoView(yaw, pitch, hfov / 1.03);
    }

    /**
     * Zoom out 3%
     */
    public void zoomOut() {
        gotoView(yaw, pitch, hfov * 1.03);
    }

    /**
     * Tilt up 5 degrees
     */
    public void panUp() {
        gotoView(yaw, pitch + hfov / pan_steps, hfov);
    }

    /**
     * Tilt down 5 degrees
     */
    public void panDown() {
        gotoView(yaw, pitch - hfov / pan_steps, hfov);
    }

    /**
     * Pan left 5 degrees
     */
    public void panLeft() {
        gotoView(yaw - hfov / pan_steps, pitch, hfov);
    }

    /**
     * Pan right 5 degrees
     */
    public void panRight() {
        gotoView(yaw + hfov / pan_steps, pitch, hfov);
    }

}
