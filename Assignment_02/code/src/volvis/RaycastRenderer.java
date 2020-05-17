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
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;


/**
 *
 * @author michel
 *  Edit by Anna Vilanova & Nicola Pezzotti
 */

// This is a very important class where you have to implement most of your work

public class RaycastRenderer extends Renderer implements TFChangeListener {

	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//in this function we update the "image" attribute using the slicing technique
    void slicer(double[] viewMatrix) {
	    // we start by clearing the image
	    resetImage();
	
	    // vector uVec and vVec define the view plane, 
	    // perpendicular to the view vector viewVec which is going from the view point towards the object
	    double[] viewVec = new double[3];
	    double[] uVec = new double[3];
	    double[] vVec = new double[3];
	    getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);
	  
	    // compute the volume center
	    double[] volumeCenter = new double[3];
	    computeVolumeCenter(volumeCenter);
	    
	    // Here will be stored the 3D coordinates of every pixel in the plane 
	    double[] pixelCoord = new double[3];
	   	
	    // sample on a plane through the origin of the volume data
	    double max = volume.getMaximum();
	    TFColor voxelColor = new TFColor();
            TFColor colorAux;
	    //Iterate on every pixel
	    for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	        	//pixelCoord now contains the 3D coordinates for pixel (i,j)
	        	computePixelCoordinates(pixelCoord,volumeCenter,uVec,vVec,i,j);
		            
	            //pixelCoord now contains the 3D coordinates of the pixels (i,j)
	            //we now have to get the value for the in the 3D volume for the pixel
	            //we can use a nearest neighbor implementation like this:
	            int val = volume.getVoxelNN(pixelCoord);

	            		
	            //you have to implement the function getVoxelLinearInterpolated in Volume.java
	            //in order to complete the assignment
	            //int val = volume.getVoxelLinearInterpolate(pixelCoord); //and then use this line
	            
	            
	            // Map the intensity to a grey value by linear scaling
	            voxelColor.r = val/max;
	            voxelColor.g = voxelColor.r;
	            voxelColor.b = voxelColor.r;
	            
	            // the following instruction makes intensity 0 completely transparent and the rest opaque
	               voxelColor.a = val > 0 ? 1.0 : 0.0;   
                       
	            // Alternatively, apply the transfer function to obtain a color using the tFunc attribute
	            // colorAux= tFunc.getColor(val);
                    //voxelColor.r=colorAux.r;voxelColor.g=colorAux.g;voxelColor.b=colorAux.b;voxelColor.a=colorAux.a; 
                    // You can also simply use voxelColor = tFunc.getColor(val); However then you copy by reference and this means that if you change 
                    // voxelColor you will be actually changing the transfer function
	            
	            //BufferedImage expects a pixel color packed as ARGB in an int
	            //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
	            int pixelColor = computeImageColor(voxelColor.r,voxelColor.g,voxelColor.b,voxelColor.a);
	            image.setRGB(i, j, pixelColor);
	        }
	    }
	}
    

  

    	//traceRay is called for every pixel of the image. Based on the number of samples, it finds the vector where the highest intensity occurs. 
       //We are dividing this intensity,which is returned by the interpollation function,with 255. as this is the highest value of rgb.Based on that value the color pixel is returned as an integer. 
    int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {

        //compute the increment and the number of samples. 
        //The samples vector indicates the coordinates of the next sample point we are going to consider in our computations, in the direction of the ray vector.
        //Ray vector points out the direction to which the viewer is looking at the object.
        double[] samples = new double[3];
        VectorMath.setVector(samples, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        //number of samples is computed as the whole 'length' of the object divided by the sample step.
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
        
        double maximum = 0;
        do {
            //devide the voxel value with 255 to normalize it, as 255 is the maximum rgb value in RGB scaling.
            double value = volume.getVoxelLinearInterpolate(currentPos)/255.;
            
            if (value > maximum) {
                maximum = value;
            }
            for (int i = 0; i < 3; i++) {
                currentPos[i] += samples[i];
            }
            nrSamples--;
        }while (nrSamples > 0);//while(VectorMath.distance(entryPoint, exitPoint)>VectorMath.distance(entryPoint, currentPos)) ;//

        double alpha;
        double r, g, b;
        if (maximum > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = maximum;
        int color = computeImageColor(r,g,b,alpha);

        return color;
    }
    
    //Composite mode is implemented here.Also the function is called when we are in transfer 2D mode (indirect vr)
    //Shading is also added based on Phong model when user selects this option
    int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        double[] lightVector = new double[3];
        double[] halfVector = new double[3];
        //the light vector is directed toward the view point (which is the source of the light)
        //half vector is used to speed up the phong shading computation see slides
        getLightVector(lightVector,halfVector,rayVector);
        
           

        // initializing color 
        double r=0;
        double g=0;
        double b=0;
        double alpha=0; 
        
        
        //compute increments along the ray for each step in a back to front mode
        double[] increments = new double[3];
        
        VectorMath.setVector(increments, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);
        //compute the number of samples that will be computed. the first one is the exit point!(back to front mode)
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //Using  back-to-front composition, therefore the current position is initialized as the exit point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);
       
        
       //instantiate Transfer function color object 
        TFColor voxel_color = new TFColor();
        TFColor colorAux = new TFColor();
        do {
        	//Gets the value  in the current position by calling the tri-linear interpolation function in volume class
            int value = volume.getVoxelLinearInterpolate(currentPos);
            //computes the gradient of this position
            VoxelGradient gradient = gradients.getGradient(currentPos);
            //Updates the color and the opacity based on the current mode selection
            if (compositingMode) {
                colorAux = tFunc.getColor(value);
                voxel_color.r =colorAux.r;
                voxel_color.g =colorAux.g;
                voxel_color.b =colorAux.b;
                //opacity
                voxel_color.a =colorAux.a;
            }       
            
 
            if (tf2dMode) {
                colorAux = tFunc2D.color;
                voxel_color.r =colorAux.r;
                voxel_color.g =colorAux.g;
                voxel_color.b =colorAux.b;
                voxel_color.a =colorAux.a;
                 voxel_color.a=tFunc2D.color.a;
                //computes the sample point opacity based on Levoy equation.
               voxel_color.a *= computeLevoyOpacity(tFunc2D.baseIntensity, tFunc2D.radius, value, gradient.mag);
            }
            //shading mode has to be applied after selecting composite mode or tf2dMode.DO NOT MOVE this if statement.
             if (shadingMode) {
                if (voxel_color.a > 0.0) {
                    //based on Phong shading function a color is assigned to the sample point
                    colorAux= computePhongShading(voxel_color, gradient, lightVector, halfVector);
                    voxel_color.r =colorAux.r;
                    voxel_color.g =colorAux.g;
                    voxel_color.b =colorAux.b;
                    voxel_color.a =colorAux.a;
                }
            }
            // Compute the composition with the back-to-front algorithm: C_i'=C_i+(1-A_i)C_(i-1)'
            r = voxel_color.a * voxel_color.r + (1.0 - voxel_color.a) * r;
            g = voxel_color.a * voxel_color.g + (1.0 - voxel_color.a) * g;
            b = voxel_color.a * voxel_color.b + (1.0 - voxel_color.a) * b;
            alpha = voxel_color.a + (1.0 - voxel_color.a) * alpha;
            

            //update the current position
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);
        
        
        //computes the color
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
    //computes the Opacity of the point based on the Levoy function. It receives as input, the material value f_v, the voxel value f_x , the radius r and the gradient magnitude.
    //The above are instanciated in the TransferFunction2D object and can be modified by the user.
       public double computeLevoyOpacity(double f_v, double r, double f_x, double gradMagnitude_fx) {

        double opacity_x;
        
        if (gradMagnitude_fx == 0.0 && f_v == f_x) {
            opacity_x = 1.0;
        } 
        else if (gradMagnitude_fx > 0.0 && (f_x - r * gradMagnitude_fx <= f_v) && (f_v <= f_x + r * gradMagnitude_fx)) {
            
            opacity_x = 1.0 - Math.abs((f_v - f_x) / (Math.abs(gradMagnitude_fx) * r));
        }
        else {
            return 0;
        }

        return opacity_x;
    }
    
   
       
       //the function is used to add shading based on Phong model, which is based on the Ambient the diffuse and the Specular coefficient
       //light vector (L) and half vector (H) are inputs, that indicate the direction in which the light falls to our object (light vector)
       // and the midpoint between the viewpoint vector and the light vector
       private TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
            double[] halfVector) {

        double k_d = 0.7;
        double k_a = 0.1;
        double k_s = 0.2;
        double n_power = 10;

        double[] NormalizedGrad = new double[3];
        VectorMath.setVector(NormalizedGrad, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);
        // inner product of light and grad vector(normal vector).
        double LNProduct = VectorMath.dotproduct(NormalizedGrad, lightVector);
        
        TFColor color = new TFColor(voxel_color.r, voxel_color.g, voxel_color.b, voxel_color.a);
        double[] V=new double[3];
        double []R=new double[3];
        
         //Viewpoint vector is :V=H-L
        //R=(2N*L)n-L
        for(int i=0;i<3;i++){
            V[i]=halfVector[i]-lightVector[i];
            R[i]=2*LNProduct*NormalizedGrad[i]-lightVector[i];
        }
      
        
        if (LNProduct > 0) {
            //L_d*k_d*c_d*(L*N)
            color.r = voxel_color.r * LNProduct * k_d + k_a*voxel_color.r;
            color.g = voxel_color.g * LNProduct *k_d + k_a*voxel_color.g;
            color.b = voxel_color.b * LNProduct * k_d + k_a*voxel_color.b;
        }
        // specular = VectorMath.dotproduct(NormalizedGrad, halfVector);
        double VRProduct=VectorMath.dotproduct(V, R);
        //adding specular coefficient
        if (VRProduct > 0) {
            
            //ALTERNATIVE to compute the specular coefficient without expotential computational cost. 
            double t=VRProduct;
            
           color.r += k_s * (t/(n_power-n_power*t+t));
            color.g += k_s *(t/(n_power-n_power*t+t));
            color.b += k_s * (t/(n_power-n_power*t+t));
            
            
            //TRADITIONAL method
            //color.r += k_s * Math.pow(VRProduct, n_power);
            //color.g += k_s * Math.pow(VRProduct, n_power);
            //color.b += k_s * Math.pow(VRProduct, n_power);
        }
        color.r = color.r > 1.0 ? 1.0 : color.r;
        color.g = color.g > 1.0 ? 1.0 : color.g;
        color.b = color.b > 1.0 ? 1.0 : color.b;
        
        return color;
    }
       
       // As an extension we could use Blinn Shading method instead of traditional Phong. Currently the function
       //is never called.
       private TFColor computeBlinnShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
            double[] halfVector) {

        double k_d = 0.7;
        double k_a = 0.1;
        double k_s = 0.2;
        double n_power = 10;

        double[] NormalizedGrad = new double[3];
        VectorMath.setVector(NormalizedGrad, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);
        // inner product of light and grad vector(normal vector).
        double LNProduct = VectorMath.dotproduct(NormalizedGrad, lightVector);
        
        TFColor color = new TFColor(voxel_color.r, voxel_color.g, voxel_color.b, voxel_color.a);
        double[] V=new double[3];
        double []R=new double[3];
        
         //Viewpoint vector is :V=H-L
        //R=(2N*L)n-L
        for(int i=0;i<3;i++){
            V[i]=halfVector[i]-lightVector[i];
            R[i]=2*LNProduct*NormalizedGrad[i]-lightVector[i];
        }
      
        
        if (LNProduct > 0) {
            //L_d*k_d*c_d*(L*N)
            color.r = voxel_color.r * LNProduct * k_d + k_a*voxel_color.r;
            color.g = voxel_color.g * LNProduct *k_d + k_a*voxel_color.g;
            color.b = voxel_color.b * LNProduct * k_d + k_a*voxel_color.b;
        }
        double specular = VectorMath.dotproduct(NormalizedGrad, halfVector);
        
        //adding specular coefficient
        if (specular > 0) {
                        
            //TRADITIONAL method
            color.r += k_s * Math.pow(specular, n_power);
            color.g += k_s * Math.pow(specular, n_power);
            color.b += k_s * Math.pow(specular, n_power);
        }
        color.r = color.r > 1.0 ? 1.0 : color.r;
        color.g = color.g > 1.0 ? 1.0 : color.g;
        color.b = color.b > 1.0 ? 1.0 : color.b;
        
        return color;
    }   
    
    void raycast(double[] viewMatrix) {

    	//data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];

        // ray parameters
         int increment;
        double sampleStep;
       
        
        if (interactiveMode) {
            increment = 2;
            sampleStep = 2.0;
            
        } 
        else{
            increment = 1;
            sampleStep = 1.0;
        }
      
        // reset the image to black
        resetImage();
        // compute the view plane and the view vector that is used to compute the entry and exit point of the ray the viewVector is pointing towards the camera
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);
        
        //The ray is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0]=-viewVec[0];rayVector[1]=-viewVec[1];rayVector[2]=-viewVec[2];
        
        // We use orthographic projection. Viewer is far away at the infinite, all pixels have the same rayVector.
        
        // ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j += increment) {
            for (int i = 0; i < image.getWidth(); i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	computePixelCoordinatesBehind(pixelCoord,viewVec,uVec,vVec,i,j);
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    if (compositingMode || tf2dMode) {
                        val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (mipMode) {
                        val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }
		
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunction2D tFunc2D;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;

    //Do NOT modify this function
    int computeImageColor(double r, double g, double b, double a){
		int c_alpha = 	a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = 	r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = 	g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = 	b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = getColorInteger(c_red,c_green,c_blue,c_alpha);
        return pixelColor;
	}
    //Do NOT modify this function    
    public void resetImage(){
    	for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	            image.setRGB(i, j, 0);
	        }
	    }
    }
    //Do NOT modify this function
    void getLightVector(double[] lightVector, double[] halfVector, double[] viewVec){
    	VectorMath.setVector(lightVector, viewVec[0], viewVec[1], viewVec[2]);
    	for (int i=0; i<3; i++) {
            halfVector[i] = viewVec[i] + lightVector[i];
        }
        double l = VectorMath.length(halfVector);
        for (int i=0; i<3; i++) {
            halfVector[i] /= l;
        }
    }
   
    //used by the slicer
    //Do NOT modify this function
    void getViewPlaneVectors(double[] viewMatrix, double viewVec[], double uVec[], double vVec[]) {
            VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
	    VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
	    VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
	}
    
    //used by the slicer	
    //Do NOT modify this function
	void computeVolumeCenter(double volumeCenter[]) {
		VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
	}
	
    //used by the slicer
    //Do NOT modify this function
	void computePixelCoordinates(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], int i, int j) {
            // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
		int imageCenter = image.getWidth()/2;
		pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
	}
    //Do NOT modify this function
	void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
		int imageCenter = image.getWidth()/2;
                // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
                // to the diaganal to make sure I am far away enough.
                double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
		pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
	}
	
    //Do NOT modify this function
    public int getColorInteger(int c_red, int c_green, int c_blue, int c_alpha) {
    	int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    	return pixelColor;
    } 
    //Do NOT modify this function
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }
    //Do NOT modify this function
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
       
        // Initialize transferfunction 
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tFunc2D= new TransferFunction2D((short) (volume.getMaximum() / 2), 0.2);
        tfEditor2D = new TransferFunction2DEditor(tFunc2D,volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }
    //Do NOT modify this function
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    //Do NOT modify this function
    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    //Do NOT modify this function
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    //Do NOT modify this function
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    public void setQuickMode(boolean mode) {
        interactiveMode = mode;
        changed();
    }
    //Do NOT modify this function
    public void setMIPMode() {
        setMode(false, true, false, false);
    }
    //Do NOT modify this function
    public void setSlicerMode() {
        setMode(true, false, false, false);
    }
    //Do NOT modify this function
    public void setCompositingMode() {
        setMode(false, false, true, false);
    }
    //Do NOT modify this function
    public void setTF2DMode() {
        setMode(false, false, false, true);
    }
    //Do NOT modify this function
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        changed();
    }
    //Do NOT modify this function
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }
    //Do NOT modify this function
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    //Do NOT modify this function
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) < 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
 
  
   
    
    //Do NOT modify this function
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    //Do NOT modify this function
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
    //Do NOT modify this function
    @Override
    public void visualize(GL2 gl) {

        double[] viewMatrix = new double[4 * 4];
        
        if (volume == null) {
            return;
        }
        	
         drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
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

    //Do NOT modify this function
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
