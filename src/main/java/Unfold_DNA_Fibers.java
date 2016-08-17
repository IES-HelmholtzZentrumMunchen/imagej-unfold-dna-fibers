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
import java.util.Vector;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
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
	public int width = 5;

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
			Vector<ImagePlus> fibers = this.process();
			
			// TODO display the output
		}
	}

	/**
	 * Show the dialog box for input parameters.
	 * @return True if the dialog box has been filled and accepted, false otherwise.
	 */
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("DNA Fibers - extract and unfold");
	
		gd.addNumericField("width", this.width, 0);
	
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
	
		this.width = (int)gd.getNextNumber();

		return true;
	}
	
	private boolean showAndCheckDialog() {
		// Call dialog
		boolean notCanceled = this.showDialog();
		
		// Check parameters
		while(notCanceled && !this.checkWidth(this.width)) {
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
	public Vector<ImagePlus> process() {
		// TODO process
		return new Vector<ImagePlus>();
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
		
	}
}
