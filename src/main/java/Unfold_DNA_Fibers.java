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

import java.awt.Button;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.Arrays;
import java.util.Vector;

import ij.CompositeImage;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.FloatPolygon;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.ContrastEnhancer;
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
	protected ImagePlus image = null;
	
	/** The input ROIs */
	protected Vector<Roi> rois = null;

	/** The width of the fibers */
	protected int radius = 4;
	
	/** Group the unfolded fibers on one image output instead of many outputs (default: true). */
	protected boolean groupFibers = true;
	
	/** X-axis margin for output fibers. */
	protected int groupMargin = 5;
	
	/** Additional space for text label on grouped output. */
	protected int groupLabelSpace = 20;
	
	/** Space between fibers on grouped output. */
	protected int groupFiberSpace = 5;
	
	/** Output the files in the specified folder instead of in current ImageJ session (default: false). */
	protected boolean output = false;
	
	/** Output path used to save files when user specified it. */
	protected String outputPath = "";

	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		// Get inputs
		this.image = imp;
		this.rois = new Vector<Roi>();
		
		// Gather the ROIs either in the manager or on the image
		RoiManager manager = RoiManager.getInstance();
		if (manager == null) {
			if (this.image.getRoi() != null)
				this.rois.add(this.image.getRoi());
		}
		else
			this.rois.addAll(Arrays.asList(manager.getRoisAsArray()));
		
		// Check if there are actually ROIs
		if (this.rois.size() == 0) {
			IJ.showMessage("This command requires at least one ROI!");
			return DONE;
		}
		
		// Check validity of ROIs
		if (!this.checkRois(this.rois)) {
			IJ.showMessage("Only straight, segmented and freehand lines are accepted as ROIs!");
			return DONE;
		}
		
		// Finish the setup
		return DOES_8G | DOES_16 | DOES_32 | NO_CHANGES;
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
	private boolean checkRois(Vector<Roi> rois) {
		for (Roi roi : rois) {
			if (roi != null && 
			roi.getType() != Roi.LINE && 
			roi.getType() != Roi.POLYLINE &&
			roi.getType() != Roi.FREELINE) {
				return false;
			}
		}
		
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
	 * Check the validity of the output path.
	 * @param outputPath Parameter value to check.
	 * @return True if the output path is valid, false otherwise.
	 */
	private boolean checkOutputPath(String outputPath) {
		File path = new File(outputPath);

		if (path.exists() && path.isDirectory())
			return true;
		else
			return false;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
		if (this.showAndCheckDialog()) {
			Vector<UnfoldedFiber> fibers = this.process();

			// Create output group image (if asked)
			ImagePlus groupOutput = null;
			
			if (this.groupFibers) {
				int maxPixelsLength = this.getMaximalFiberLengthInPixels(fibers);
				
				// To create group output, we need the height (2 x nb of fibers x (2xradius+1)) 
				// and the width (maximal width of individual fibers).
				groupOutput = IJ.createHyperStack(
						"Fibers of "+this.image.getTitle(), 
						maxPixelsLength+2*this.groupMargin+this.groupLabelSpace, 
						fibers.size()*(2*this.radius+1+2*this.groupFiberSpace), this.image.getNChannels(),
						this.image.getNSlices(), this.image.getNFrames(), this.image.getBitDepth());
				
				if (!this.output)
					groupOutput.show();
			}
			
			// If composite image, try to get the channel's color
			CompositeImage composite = null;
			if (this.image.isComposite())
				composite = (CompositeImage)this.image;
			
			// Get image and plot from each fiber
			for (int i = 0; i < fibers.size(); i++) {
				// Display fiber image or group output image
				if (this.groupFibers)
					this.copyFiberIntoGroup(groupOutput, fibers.get(i), i);
				else { // !this.groupFibers
					if (this.output)
						IJ.save(fibers.get(i).fiberImage, this.outputPath+fibers.get(i).fiberImage.getTitle()+".zip");
					else // !this.output
						fibers.get(i).fiberImage.show();
				}
				
				// Display the profiles with the plot GUI
				Plot plot = new Plot("Profiles #"+IJ.d2s(i+1,0), "Length ["+this.image.getCalibration().getXUnit()+"]", "Intensity level [a.u.]");
				
				for (int c = 0; c < fibers.get(i).fiberProfiles.size(); c++) {
					if (this.image.isComposite()) { // Match channel's color with profile's color
						composite.setC(c+1);
						plot.setColor(composite.getChannelColor());
					}
					
					plot.addPoints(fibers.get(i).profilesAbscissa, fibers.get(i).fiberProfiles.get(c), Plot.LINE);
					plot.draw();
				}
				
				if (this.output) {
					IJ.save(plot.getImagePlus(), this.outputPath+plot.getTitle()+".png");
					plot.getResultsTable().save(this.outputPath+plot.getTitle()+".csv");
				}
				else // !this.output
					plot.show();
			}
			
			if (this.groupFibers && this.output)
				IJ.save(groupOutput, this.outputPath+groupOutput.getTitle()+".zip");
		}
	}
	
	/**
	 * Copy the fiber image into group image.
	 * @param groupImage Target group image. 
	 * @param fiber Fiber to copy.
	 * @param index Index of fiber to copy.
	 */
	private void copyFiberIntoGroup(ImagePlus groupImage, UnfoldedFiber fiber, int index) {
		int y0 = index*(2*this.radius+1+2*this.groupFiberSpace)+this.groupFiberSpace;
		
		for (int c = 1; c <= groupImage.getNChannels(); c++) {
			groupImage.setC(c);
			fiber.fiberImage.setC(c);
			
			ImageProcessor gip = groupImage.getChannelProcessor();
			ImageProcessor fip = fiber.fiberImage.getChannelProcessor();
			
			for (int y = 0; y < fip.getHeight(); y++) {
				for (int x = 0; x < fip.getWidth(); x++) {
					gip.setf(x+this.groupMargin+this.groupLabelSpace, y0+y, fip.getf(x, y));
				}
			}
			
			gip.setValue((int)Math.pow(2, 8)-1);
			gip.drawString("#"+IJ.d2s(index+1, 0), 0, y0+this.radius+8);
			
		}
	}
	
	/**
	 * Get the maximal number of pixels of fiber images.
	 * @param fibers Unfolded fibers.
	 * @return The maximal length in pixels (as integer).
	 */
	private int getMaximalFiberLengthInPixels(Vector<UnfoldedFiber> fibers) {
		int maxPixelsLength = 0;
		
		for (UnfoldedFiber fiber : fibers) {
			if (fiber.fiberImage.getWidth() > maxPixelsLength)
				maxPixelsLength = fiber.fiberImage.getWidth();
		}
		
		return maxPixelsLength;
	}

	/**
	 * Show the dialog box for input parameters.
	 * @return True if the dialog box has been filled and accepted, false otherwise.
	 */
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("DNA Fibers - extract and unfold");

		// Add output selection
		gd.addStringField("output path", this.outputPath, 30);
		Button firstButton = new Button("Select folder");
        firstButton.addActionListener(new ActionListener() {
            // Execute when button is pressed
            public void actionPerformed(ActionEvent e) {
                // Display a browser
                String path = IJ.getDirectory("Select output path");

                // Display file path in text field
                if (path != null) {
                    Vector<TextField> stringFields = gd.getStringFields();
                    stringFields.get(0).setText(path);
                }
            }
        });
        gd.add(firstButton);
        
        // Add other components
		gd.addNumericField("radius", this.radius, 0, 5, "pixels");
		gd.addCheckbox("group fibers", this.groupFibers);
		gd.addCheckbox("auto save files to output path", this.output);
	
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
	
		this.outputPath  = gd.getNextString();
		this.radius      = (int)gd.getNextNumber();
		this.groupFibers = gd.getNextBoolean();
		this.output      = gd.getNextBoolean();

		return true;
	}
	
	private boolean showAndCheckDialog() {
		// Call dialog
		boolean notCanceled = this.showDialog();
		
		// Check parameters
		while ( notCanceled &&
				(!this.checkWidth(this.radius) ||
				(this.output && !this.checkOutputPath(this.outputPath))) ) {
			IJ.showMessage("Width must be strictly positive!\nOutput path must be a valid and existing path to a folder!");
			notCanceled = this.showDialog();
		}
		
		return notCanceled;
	}

	/**
	 * Extract and unfold DNA fibers selected by ROIs.
	 * 
	 * For now, this method support only one ROI.
	 * 
	 * For unfolding, the tangents of the median lines of the fibers
	 * need to be computed first. Then multi-channel images are sampled
	 * along the normal of the tangents, with a defined radius. New image
	 * containing the unfolded fiber and profiles are created from these
	 * samples.
	 *
	 * @return A vector of images containing extracted and unfolded DNA fibers.
	 */
	public Vector<UnfoldedFiber> process() {
		Vector<UnfoldedFiber> fibers = new Vector<UnfoldedFiber>();
		
		// Each ROI is processed separately
		for (int i = 0; i < this.rois.size(); i++) {
			Vector<Point2D[]> normals = this.computeNormals(this.rois.get(i));
	
			// Sample intensity along normals (2*radius + 1) and compute
			// the profile (take the maximal value of each column).
			// Put the intensities into a new image of unfolded fiber
			// and the maximal intensity in a table.
			
			ImagePlus unfoldedFiber = IJ.createHyperStack( // Initialize the image of unfolded fiber
					"Fiber #"+IJ.d2s(i+1,0), normals.size(), 2*this.radius+1, this.image.getNChannels(),
					this.image.getNSlices(), this.image.getNFrames(), this.image.getBitDepth());
			
			Vector<double[]> profiles = new Vector<double[]>();
			double[] profilesAbscissa = new double[normals.size()];
			
			
			// Go through each channels
			for (int c = 1; c <= this.image.getNChannels(); c++) {
				this.image.setC(c);
				ImageProcessor ip = this.image.getChannelProcessor();
				
				unfoldedFiber.setC(c);
				ImageProcessor fp = unfoldedFiber.getChannelProcessor();
				
				ip.setInterpolate(true);
				ip.setInterpolationMethod(ImageProcessor.BICUBIC);
				
				double[] profile = new double[normals.size()];
				
				// Go through each point of ROI
				for (int x = 0; x < normals.size(); x++) {
					profilesAbscissa[x] = x * this.image.getCalibration().pixelWidth;
					
					Point2D[] normal = normals.get(x);
					profile[x] = 0.;
					
					for (int s = -this.radius, y = 0; s <= this.radius; s++, y++) {
						double interpolatedValue = ip.getInterpolatedPixel( // Interpolate pixel value
								normal[0].getX() + s*normal[1].getX(), 
								normal[0].getY() + s*normal[1].getY());
	
						fp.setf(x,y,(float)interpolatedValue); // Fill image of unfolded fiber
						
						if (Double.compare(profile[x], interpolatedValue) < 0) // Get maximal value for profile
							profile[x] = interpolatedValue;
					}
				}
				
				profiles.add(profile);
			}
			
	
			// If image is composite, copy the channel's color scheme
			if (this.image.isComposite()) {
				unfoldedFiber.setDisplayMode(IJ.COMPOSITE);
				CompositeImage composite = (CompositeImage)unfoldedFiber;
				composite.copyLuts(this.image);
				composite.setMode(CompositeImage.COLOR); composite.setMode(CompositeImage.COMPOSITE); // This trick is needed to force the display
			}
			
			
			fibers.add(new UnfoldedFiber(unfoldedFiber, profiles, profilesAbscissa));
		}
		
		return fibers;
	}
	
	/**
	 * Compute the normals from the points of a ROI.
	 * 
	 * The ROI is supposed to be a line, a polyline or a freeline.
	 * The output is a vector of arrays, containing two elements:
	 * the first one is the current point and the second one is 
	 * the associated normal.
	 * 
	 * @param roi The input ROI.
	 * @return The points in image and associated normals of the given ROI.
	 */
	protected Vector<Point2D[]> computeNormals(Roi roi) {
		Vector<Point2D[]> normals = new Vector<Point2D[]>();
		
		/*RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();*/
		
		FloatPolygon points = roi.getInterpolatedPolygon();
		
		// Computing the tangent needs 5 points; therefore 2 points at the beginning and
		// 2 points at the end will not be used for the unfolding. This issue can be
		// overcome by extrapolation (e.g. linear).
		// In our case, it represents only a fiber path shortened by 2 pixels at the
		// beginning and 2 pixels at the end. So we omit a correction for that issue.
		
		// Go through each point
		for (int i = 2; i < points.npoints-2; i++) {
			// Instead using finite differences, the tangent vector is computed at 
			// point p with a least-square linear fit 5 points in total (the central
			// point 2 points before and 2 points after) in order to avoid numerical
			// issues.
			double x1 = points.xpoints[i-2], y1 = points.ypoints[i-2];
			double x2 = points.xpoints[i-1], y2 = points.ypoints[i-1];
			double x3 = points.xpoints[i],   y3 = points.ypoints[i];
			double x4 = points.xpoints[i+1], y4 = points.ypoints[i+1];
			double x5 = points.xpoints[i+2], y5 = points.ypoints[i+2];
			
			double slope = ( x1*y2 - 4.*x1*y1 + x2*y1 + x1*y3 - 4.*x2*y2 + x3*y1 + x1*y4 + x2*y3 + x3*y2 + x4*y1 + 
							 x1*y5 + x2*y4 - 4.*x3*y3 + x4*y2 + x5*y1 + x2*y5 + x3*y4 + x4*y3 + x5*y2 + x3*y5 - 
							 4.*x4*y4 + x5*y3 + x4*y5 + x5*y4 - 4.*x5*y5 ) / 
						   ( 2. * ( - 2.*x1*x1 + x1*x2 + x1*x3 + x1*x4 + x1*x5 - 2.*x2*x2 + x2*x3 + x2*x4 + x2*x5 - 
								   2.*x3*x3 + x3*x4 + x3*x5 - 2.*x4*x4 + x4*x5 - 2.*x5*x5 ) );
			
			Point2D normal = new Point2D.Double(1,0); // In NaN case, it means horizontal normal (vertical tangent), so initialize to (1,0)
			
			if (!Double.isNaN(slope)) {
				double x = Math.sqrt(1./(1+slope*slope)); // Parametric formula that makes unit vector
				normal = new Point2D.Double(-slope*x, x); // In 2D, orthogonal vector is unique, so closed-form solution
			}
			
			// Fix heterogeneous orientation of normal vector in order to get 
			// consistent results (e.g. image unfolding).
			if (i > 2) {
				Point2D lastNormal = normals.lastElement()[1];
				
				// Since vectors are unit, if the dot product is negative, they have opposite orientations
				if (Double.compare(normal.getX()*lastNormal.getX() + normal.getY()*lastNormal.getY(),0) < 0) {
					normal.setLocation(-normal.getX(), -normal.getY());
				}
			}
			
			/*manager.addRoi(new PointRoi(x3,y3));
			manager.addRoi(new Line(
					x3-this.radius*normal.getX(), y3-this.radius*normal.getY(),
					x3+this.radius*normal.getX(), y3+this.radius*normal.getY()));
			manager.addRoi(new PointRoi(x3+this.radius*normal.getX(), y3+this.radius*normal.getY()));*/
			
			normals.add(new Point2D[]{new Point2D.Double(x3,y3), normal});
		}
		
		// TODO smooth the slope with a moving average?
		
		return normals;
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
		public ImagePlus fiberImage = null;
		
		/**
		 * Intensity profiles of the unfolded fiber.
		 */
		public Vector<double[]> fiberProfiles = null;
		
		/**
		 * Abscissa of intensity profiles.
		 */
		public double[] profilesAbscissa = null;
		
		/**
		 * Default constructor
		 */
		UnfoldedFiber() {
			
		}
		
		/**
		 * Constructor
		 * @param fiberImage
		 * @param fiberProfiles
		 */
		UnfoldedFiber(ImagePlus fiberImage, Vector<double[]> fiberProfiles, double[] profilesAbscissa) {
			this.fiberImage = fiberImage;
			this.fiberProfiles = fiberProfiles;
			this.profilesAbscissa = profilesAbscissa;
		}
	}
}
