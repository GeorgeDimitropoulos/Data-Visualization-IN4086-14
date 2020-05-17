/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 *  Modified by Anna Vilanova
 */
public class GradientVolume {

	
	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	//Compute the gradient of contained in the volume attribute and save it into the data attribute
	//
	//This is a lengthy computation and is performed only once (have a look at the constructor GradientVolume) 

    
    //this function calculates the gradients of the vertices of the cubes. 
    //We need those gradients, as we will use them to find the normal vector of the 
    //isosurfaces(triangles) through interpolation.
    private void compute() {

        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
          for (int z=1; z<dimZ-1; z++) {
            for (int y=1; y<dimY-1; y++) {
                for (int x=1; x<dimX-1; x++) {
                    float gx = (volume.getVoxel(x+1, y, z) - volume.getVoxel(x-1, y, z))/2.0f;
                    float gy = (volume.getVoxel(x, y+1, z) - volume.getVoxel(x, y-1, z))/2.0f;
                    float gz = (volume.getVoxel(x, y, z+1) - volume.getVoxel(x, y, z-1))/2.0f;
                    
                    //creating a voxelgradient instance.
                    VoxelGradient grad = new VoxelGradient(gx, gy, gz);
                    //updating data[] with the gradient computed
                    setGradient(x, y, z, grad);
                }
            }
        }
    }
    	
    //This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    //the result is given at result. You can use it to tri-linearly interpolate the gradient
    private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {
        //f(s)=(D-d)f(x1)+d*f(x2) , where D is the distance between x1 and x2 and d is the distance between x1 and s.
        result.x = (1.0f - factor)*g0.x + factor*g1.x;
        result.y = (1.0f - factor)*g0.y + factor*g1.y;
        result.z = (1.0f - factor)*g0.z + factor*g1.z;
        
        //magnitude
        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
    
    }
    
 
    // This function returns the gradient at position coord, using trilinear interpolation between the gradients of the cube.   
    public VoxelGradient getGradient(double[] coord) {

        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        //distance between coord and its nearest cube vertex in x,y,z axis
        float d_x = (float) coord[0] - x;
        float d_y = (float) coord[1] - y;
        float d_z = (float) coord[2] - z;

        
        
        //interpolation between gradients in (x,y,z) and (x+1,y,z)
        VoxelGradient c_00 = new VoxelGradient();
        interpolate(getGradient(x, y, z), getGradient(x+1, y, z), d_x, c_00);
        
        //interpolation between gradients in (x,y+1,z) and (x+1,y+1,z)
        VoxelGradient c_10 = new VoxelGradient();
        interpolate(getGradient(x, y+1, z), getGradient(x+1, y+1, z), d_x, c_10);
        
        //interpolation between gradients in (x,y,z+1) and (x+1,y,z+1)
        VoxelGradient c_01 = new VoxelGradient();
        interpolate(getGradient(x, y, z+1), getGradient(x+1, y, z+1), d_x, c_01);
        
        //interpolation between gradients in (x,y+1,z+1) and (x+1,y+1,z+1)
        VoxelGradient c_11 = new VoxelGradient();
        interpolate(getGradient(x, y+1, z+1), getGradient(x+1, y+1, z+1), d_x, c_11);

        //interpolation between the gradients already computed in the midpoints c_00 and c_10 in the y axis
        VoxelGradient c_0 = new VoxelGradient();
        interpolate(c_00, c_10, d_y, c_0);
 
        //interpolation between the gradients already computed in the midpoints c_00 and c_10 in the y axis        
        VoxelGradient c_1 = new VoxelGradient();
        interpolate(c_01, c_11, d_y, c_1);
        
        //we finally compute the gradient in coord by interpolating the computed gradients in c_0 and c_1.
        VoxelGradient c = new VoxelGradient();
        interpolate(c_0, c_1, d_z, c);
        return c;

       
    }
    
    
    
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	
    //Returns the maximum gradient magnitude
    //
    //The data array contains all the gradients, in this function you have to return the maximum magnitude of the vectors in data[] 
 
    //Do NOT modify this function
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    	
	
	
	//Do NOT modify this function
	public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

	//Do NOT modify this function
	public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

  
	//Do NOT modify this function: Basically calculates the Nearest Neighbor interpolation for the gradient
    public VoxelGradient getGradient2(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        return getGradient(x, y, z);
    }

	//Do NOT modify this function
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

	//Do NOT modify this function
    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }
    
	//Do NOT modify this function
    public VoxelGradient getVoxel(int i) {
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

	//Do NOT modify this attributes
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
    
    //If needed add new attributes here:
}
