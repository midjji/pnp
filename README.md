# pnp
A RANSAC and BA based pnp wrapper for the Lambdatwist p3p solver. See: 
http://openaccess.thecvf.com/content_ECCV_2018/html/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.html
Take a look at the benhmark at https://github.com/midjji/lambdatwist-p3p for comparisons to other methods. 


This is a simple but fast pnp solver which for most problems does not require you to think at all. 
By default it will work well in most cases, but have a look in the parameters file.

So I keep seeing people using very heavy robust loss functions for pnp, instead of a well considered ransac loop. 
In my experience the latter wins in every case except very high inlier noise, something which doesnt actually happen in practice for pnp problems. Avoid Opencvs Epnp in particular. 



Pnp estimates the pose of a camera in the coordinate system of the 3d points used. 
So if you have a pinhole camera observations y_n of 3d points x_world, the pnp(yn, x_world) will give you the pose which matches the observation model

x_c = x_camera = P_cw x_w

y_n = project(x_c0/x_c2, x_c1/x_c2)

If you have a camera with known camera intrinsics you need to compute the corresponding pinhole normalized coordinates, 

e.g. if you have the linear intrinsics K and the observation model y_rowcol (pixels)= K project(Pcw x_world)

then compute yn = inv(K)y_rowcol

if the lens model is non-linear you will need to do invert the lens model, at least for where you are interested in it.  The best way to do this if the intrinsics are calibrated pre use as is required for pnp is to offline estimate a regularized 2d spline surface over the image and use as a lookup table. 

If the camera is offset from the vehicle, i.e. you are using the observation model y_rowcol (pixels)= distort(K, th, project(Pcv Pvw x_world)) and the extrinsics (Pcv =P_camera_vehicle) are known, then normalize the coordinates as before and compute pnp Pcw, followed by finding the pose of vehicle in world as Pvw = inv(Pcv)Pcw
