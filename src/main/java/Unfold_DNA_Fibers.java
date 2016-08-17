/*
 * Manipulate and analyse DNA fibers data
 * This plugin extracts and unfold the DNA fibers selected by a curve ROI
 * Copyright (C) 2016  Julien Pontabry (Helmholtz IES)

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.awt.Point;
import java.awt.geom.Point2D;
import java.util.Vector;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Arrow;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

/**
 * Plugin for DNA fibers extraction and unfolding.
 *
 * After DNA fibers detection (manually or automatically), the fibers are
 * located by their centerline (represented by a line ROI). They can be
 * unfolded (to make them straight) with this plugin.
 *
 * @author Julien Pontabry
 */
public class Unfold_DNA_Fibers implements PlugInFilter {
	/** The input image */
	protected ImagePlus image;
	
	/** The input ROIs */
	protected Roi roi;

	/** The width of the fibers */
	public int radius = 2;

	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		// Get inputs
		this.image = imp;
		this.roi = this.image.getRoi();
		
		// Check validity of ROIs
		if (!this.checkRois(this.roi)) {
			IJ.showMessage("Only straight, segmented and freehand lines are accepted as ROIs!");
			return DONE;
		}
		
		// Finish the setup
		return DOES_8G | DOES_16 | DOES_32 | NO_CHANGES | ROI_REQUIRED;
	}
	
	/**
	 * Check the validity of the ROIs.
	 * 
	 * A valid ROI is not null and one of the following: line,
	 * polyline and freeline.
	 * 
	 * @param roi Input roi to check
	 * @return True if the ROIs are valid, false otherwise.
	 */
	private boolean checkRois(Roi roi) {
		if (this.roi != null && 
		this.roi.getType() != Roi.LINE && 
		this.roi.getType() != Roi.POLYLINE &&
		this.roi.getType() != Roi.FREELINE) {
			return false;
		}
		else // not null and one of the valid ROI
			return true;
	}
	
	/**
	 * Check the validity of the width parameter.
	 * @param width Parameter value to check.
	 * @return True if the width is valid, false otherwise.
	 */
	private boolean checkWidth(int width) {
		if (width <= 0)
			return false;
		else
			return true;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
		if (this.showAndCheckDialog()) {
			Vector<UnfoldedFiber> fibers = this.process();
			
			// TODO Display the output and plot the profiles
			for (UnfoldedFiber fiber : fibers) {
				fiber.fiberImage.show();
			}
		}
	}

	/**
	 * Show the dialog box for input parameters.
	 * @return True if the dialog box has been filled and accepted, false otherwise.
	 */
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("DNA Fibers - extract and unfold");
	
		gd.addNumericField("radius", this.radius, 0);
	
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
	
		this.radius = (int)gd.getNextNumber();

		return true;
	}
	
	private boolean showAndCheckDialog() {
		// Call dialog
		boolean notCanceled = this.showDialog();
		
		// Check parameters
		while(notCanceled && !this.checkWidth(this.radius)) {
			IJ.showMessage("Width must be strictly positive!");
			notCanceled = this.showDialog();
		}
		
		return notCanceled;
	}

	/**
	 * Extract and unfold DNA fibers selected by ROIs.
	 *
	 * @return A vector of images containing extracted and unfolded DNA fibers.
	 */
	public Vector<UnfoldedFiber> process() {
		Vector<UnfoldedFiber> fibers = new Vector<UnfoldedFiber>();
		
		RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();
		
		Point[] points = roi.getContainedPoints();
		
		// TODO Add points before and after curve to compensate for the point loss when fitting
		
		for (int i = 2; i < points.length-2; i++) {
			// Compute the tangent vector at point p
			// with a least-square line fit in order to
			// avoid numerical issues
			double x1 = points[i-2].getX(), y1 = points[i-2].getY();
			double x2 = points[i-1].getX(), y2 = points[i-1].getY();
			double x3 = points[i].getX(),   y3 = points[i].getY();
			double x4 = points[i+1].getX(), y4 = points[i+1].getY();
			double x5 = points[i+2].getX(), y5 = points[i+2].getY();
			
			double slope = ( x1*y2 - 4.*x1*y1 + x2*y1 + x1*y3 - 4.*x2*y2 + x3*y1 + x1*y4 + x2*y3 + x3*y2 + x4*y1 + 
							 x1*y5 + x2*y4 - 4.*x3*y3 + x4*y2 + x5*y1 + x2*y5 + x3*y4 + x4*y3 + x5*y2 + x3*y5 - 
							 4.*x4*y4 + x5*y3 + x4*y5 + x5*y4 - 4.*x5*y5 ) / 
						   ( 2. * ( - 2.*x1*x1 + x1*x2 + x1*x3 + x1*x4 + x1*x5 - 2.*x2*x2 + x2*x3 + x2*x4 + x2*x5 - 
								   2.*x3*x3 + x3*x4 + x3*x5 - 2.*x4*x4 + x4*x5 - 2.*x5*x5 ) );
			
			Point2D normal = new Point2D.Double(1,0); // In NaN case, it means vertical line (horizontal normal) with top-left origin
			
			if (!Double.isNaN(slope)) {
				double x = Math.sqrt(1./(1+slope*slope)); // Parametric formula that makes unit vector
				normal = new Point2D.Double(-slope*x, x); // In 2D, orthogonal vector is unique
			}
			
			if (i % 5 == 0) {
				manager.addRoi(new Line(
						x2-this.radius*normal.getX(), y2-this.radius*normal.getY(),
						x2+this.radius*normal.getX(), y2+this.radius*normal.getY()));
			}
			
			// TODO sample intensity along normal (2*radius + 1)
			
			// TODO 
		}
		
		return fibers;
	}

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Unfold_DNA_Fibers.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		ImagePlus image = IJ.openImage("http://imagej.net/images/blobs.gif");
		image.show();

		/*// run the plugin
		IJ.runPlugIn(clazz.getName(), "");*/
	}
	
	/**
	 * Structure for the output of the plugin.
	 * 
	 * This structure holds an image corresponding to the unfolded
	 * fibers and its associated intensity profiles.
	 * 
	 * @author Julien Pontabry
	 */
	public class UnfoldedFiber {
		/** Image of the unfolded fiber. */
		public ImagePlus fiberImage;
		
		/**
		 * Intensity profiles of the unfolded fiber.
		 *
		 * The first row contains the positions (1D
		 * coordinates) on the fiber. The further rows
		 * contain the intensity profiles.
		 */
		public double[][] fiberProfiles;
	}
}
