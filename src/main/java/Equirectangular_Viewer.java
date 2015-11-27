/* Equirectangular Viewer - Interactive equirectangular panorama Viewer for ImageJ/Fiji
 * Copyright (C) 2015 - Yili Zhao panovr@gmail.com
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
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GUI;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import javax.swing.SwingUtilities;

public class Equirectangular_Viewer implements PlugIn {

    private Viewer viewer;

    @Override
    public void run(String arg) {
        final ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) {
            IJ.beep();
            IJ.showStatus("No equirectangular image opened.");
            return;
        }

        final ColorProcessor processor = (ColorProcessor) imp.getProcessor();
        final int[] pixels = (int[]) (processor.getPixels());

        final int[][] pd = new int[imp.getHeight()][imp.getWidth()];

        for (int y = 0; y < imp.getHeight(); y++) {
            final int yw = y * imp.getWidth();
            for (int x = 0; x < imp.getWidth(); x++) {
                pd[y][x] = pixels[yw + x];
            }
        }

        viewer = new Viewer(pd);
        viewer.init();

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                createAndShowGUI(viewer);
                WindowManager.addWindow(viewer);
            }
        });
    }

    private void createAndShowGUI(final Viewer viewer) {
        viewer.setSize(400, 300);
        viewer.init();
        GUI.center(viewer);
        viewer.setVisible(true);
    }

    /**
     * Main method for debugging.
     *
     * For debugging, it is convenient to have a method that starts ImageJ,
     * loads an image and calls the plugin, e.g. after setting breakpoints.
     *
     * @param args unused
     */
    public static void main(String[] args) {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        Class<?> clazz = Equirectangular_Viewer.class;
        String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
        String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
        System.setProperty("plugins.dir", pluginsDir);

        // start ImageJ
        new ImageJ();

        // open image
        //ImagePlus image = IJ.openImage("D:/Develop/pano.jpg");
        ImagePlus image = IJ.openImage("http://fly.mpi-cbg.de/~saalfeld/Projects/download/panorama/theaterplatz_3400.jpg");
        image.show();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}
