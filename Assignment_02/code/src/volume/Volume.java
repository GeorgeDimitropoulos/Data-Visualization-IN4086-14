/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel
 *  Modified by Anna Vilanova
 */
public class Volume {
    
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
    //This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. You can use it to tri-linearly interpolate the values 
	private float interpolate(float g0, float g1, float factor) {
        float result=0;
        //this corresponds to linear interpollation function between two known points f[i],f[i+1] and a factor t:f'(x)=f[i](1-t)+f[i+1]t
        result=g0*(1-factor)+g1*factor;
        return result; 
    }
	
	//The function receives the coordinates x,y,z of a given point as an input.The goal of the function is to 
        //compute the value of this point, by performing tri-linear interpolation.
	public short getVoxelLinearInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        //we first get the coordinates x,y,z of the nearest neighboor point.
        int x = (int) Math.floor(coord[0]); 
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        //calculating t parameter (distance) between our sample point and its neighboors in the x,y and z axis
        float d_x = (float) coord[0] - x;
        float d_y = (float) coord[1] - y;
        float d_z = (float) coord[2] - z;

        //performing tri-linear interpolation. We consider a cube whose coordinates are the nearest neighboor known points in the 3d space.
        
        //interpolation between (x,y,z) neighbor point and (x+1,y,z)  
        float c0 = interpolate(getVoxel(x, y, z), getVoxel(x+1, y, z), d_x);
        
        //interpolation between (x,y,z+1) point and (x+1,y,z+1)  
        float c1 = interpolate(getVoxel(x, y, z+1), getVoxel(x+1, y, z+1), d_x);
         //interpolation between (x,y+1,z) point and (x+1,y+1,z)
        float c2 = interpolate(getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), d_x);
        
        //interpolation between (x,y+1,z+1) point and (x+1,y+1,z+1)
        float c3 = interpolate(getVoxel(x, y+1, z+1), getVoxel(x+1, y+1, z+1), d_x);
         //interpolation between c0 point and c2 point. The line connecting c0 and c2 is parrallel with the y axis
        float c4 = interpolate(c0, c2, d_y);
         //interpolation between c1 point and c3 point. The line connecting c3 and c1 is parrallel with the y axis
        float c5 = interpolate(c1, c3, d_y);
        //finally c6 is calculated as an interpolation between the previously calculated c4 and c5. Practically this is the value of the function for the intermediate point that results through interpolation.
        float c6 = interpolate(c4, c5, d_z);
            
        return (short)c6; 
    }
		
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//Do NOT modify this function
        // This function is an example and does a nearest neighbour interpolation
	public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]); 
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
    
        return getVoxel(x, y, z);
    }
	
	//Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
	//Do NOT modify this function
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
	//Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
    	int i = x + dimX*(y + dimY * z);
        return data[i];
    }
    
	//Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
    	int i = x + dimX*(y + dimY * z);
        data[i] = value;
    }
    
	//Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
	//Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }
    
	//Do NOT modify this function
    public int getDimX() {
        return dimX;
    }
    
	//Do NOT modify this function
    public int getDimY() {
        return dimY;
    }
    
	//Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

	//Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }
    
	//Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
	//Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }
    
	//Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
	//Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
