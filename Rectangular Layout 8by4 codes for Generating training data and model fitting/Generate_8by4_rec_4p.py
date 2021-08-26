import numpy as np
import pandas as pd

###Generate virtual training data 
##assign variables for crack coordinates, angle, stress intensity factors and T stress
x0_all,y0_all,b0_all,K_1_all,K_2_all,T_all=([] for i in range(6))

##Number of points and their range
#Parameters in sequence:Stress Intensity Factor I, II, T-stress,angle & Coordinates of crack(y&x)
#Number of points considered in each parameter = num_points
num_points=4
##Ranges for data parameters: 
#K_I=-300 to 1200; K_II=-110 to 130; T-stress=-60 to 22; angle=0 to 90; x_coordinate=-35 to 35; y_coordinate
K_I_range =np.linspace(0,1200,num_points)
K_II_range=np.linspace(-110,130,num_points)
T_range =np.linspace(-60,22,num_points)
Beta_range =np.linspace(0,2*np.pi,num_points)
Y_coord_range =np.linspace(-75,-6,num_points)
X_coord_range =np.linspace(-35,-6,num_points)

#Looping over to create data
for a in K_I_range:
    for b in K_II_range:
        for c in T_range:
            for angle in Beta_range:    
                for j in Y_coord_range:    
                    for i in X_coord_range:
                        x0_all.append(i)        
                        y0_all.append(j)
                        b0_all.append(angle)
                        K_1_all.append(a)
                        K_2_all.append(b)
                        T_all.append(c)
                      
#y_data = np.column_stack((x0_all, y0_all, b0_all))

#print (x0_all)
#print (y0_all)
#print (b0_all)
print (len(x0_all))

##Obtaing strain field for varying x_0 and y_0 generated above
count_2 = 0
for count_1 in range(0,len(x0_all)):
	#from the local cordinate system and beta, let's obtain global coordinate system
    x_0=x0_all[count_1]
    y_0=y0_all[count_1]
    beta = b0_all[count_1]
    #Values of stress intensity factors extracted from the research journals Table 5: pure mode_I loading with different crack position

    SIF_1=K_1_all[count_1]
    SIF_2=K_2_all[count_1]
    T=T_all[count_1]

    
    
    hor_arr=8 #Horizontal arrangement of sensors
    ver_arr=4 #vertical arrangement of sensors
    num_sensors=hor_arr*ver_arr #number of sensors are regular rectangular arrangement 
    x_dash = np.zeros(num_sensors)
    y_dash = np.zeros(num_sensors)

    x = np.zeros(num_sensors)
    y = np.zeros(num_sensors)
    phi = np.zeros(num_sensors)
    x_count = 0
    y_count = 0
    x_dist=10 #distance between sensors
    y_dist=10    


    #below is the snip for x_dash, y_dash coordinates w.r.t film C.S.
    #np.add(y_dash,y_dist) 
    count = 0
    for k in range(1,ver_arr+1):
        for l in range(1,hor_arr+1):
            x_dash[count]= 0+x_dist*l
            y_dash[count]= 0+y_dist*k
            count  = count+1   

    #obtaining global coordinate values
    for i in range(0,num_sensors):
        x[i] = x_dash[i]*np.cos(beta)-y_dash[i]*np.sin(beta)+x_0
        y[i] = y_dash[i]*np.sin(beta)+y_dash[i]*np.cos(beta)+y_0 


    
    #coordinates of the measuring points (r,phi)
    r = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)

    #Values of stress intensity factors extracted from the research journals Table 5: pure mode_I loading with different crack position
    nu = 0.33
    E = 72

    #with warnings.catch_warnings():
        #warnings.filterwarnings('error')
        
        
    np.seterr(divide='raise')
    try:
        stress_11 = ((SIF_1/(np.sqrt(2*np.pi*r))) * (np.cos(phi/2)) * (1-np.sin(phi/2)*np.sin((3*phi)/2))) - ((SIF_2/(np.sqrt(2*np.pi*r))) * np.sin(phi/2) * (2+np.cos(phi/2)*np.cos(3*phi/2))) + T
        stress_22 = ((SIF_1/(np.sqrt(2*np.pi*r))) * (np.cos(phi/2)) * (1+np.sin(phi/2)*np.sin((3*phi)/2))) + ((SIF_2/(np.sqrt(2*np.pi*r))) * np.sin(phi/2) * (np.cos(phi/2) * np.cos(3*phi/2)))
        stress_12 = ((SIF_1/(np.sqrt(2*np.pi*r))) * (np.sin(phi/2)) * (np.cos(phi/2) * np.cos((3*phi)/2))) + ((SIF_2/(np.sqrt(2*np.pi*r))) * np.cos(phi/2) * (1-(np.sin(phi/2)*np.sin(3*phi/2))))
        #print("stress completed")

    except FloatingPointError:
        #print(count_1)
        continue
    
        
    
    #strains of the specimen
    eps_specimen_11 = (stress_11-(nu*stress_22))/E
    eps_specimen_22 = (stress_22-(nu*stress_11))/E
    eps_specimen_12 = (1+nu)*stress_12/E
    #print("strain_completed")
    print(count_1)

    

    #combining all the individual arrays for the sensors placed and forming a single array 
    stack_arr = np.column_stack((eps_specimen_11,eps_specimen_22,eps_specimen_12))
    xgroup = stack_arr.flatten()
    ygroup = np.column_stack((x_0,y_0,np.cos(beta),np.sin(beta),SIF_1,SIF_2,T))


    if (count_2 == 0):
        final_xgroup = xgroup
        final_ygroup = ygroup
        count_2 = 1
    else:
        final_xgroup = np.vstack((final_xgroup, xgroup))
        final_ygroup = np.vstack((final_ygroup, ygroup))
        
x_data = final_xgroup
y_data = final_ygroup

print (np.shape(x_data),"shape of x_data") #16 sensors and 5 outputs for each sensor, so 80
print(np.shape(y_data),"shape of y_data") # crack location(x and y) and angle(so 3 columns)

#print (np.shape(x_data))
#print(np.shape(y_data))    

#print (x_data[:,:4])
#print (y_data)
df_x = pd.DataFrame(data=x_data)
df_y = pd.DataFrame(data=y_data)
df_x.to_csv('x_rec_32_data_4.csv')
df_y.to_csv('y_rec_32_data_4.csv')

