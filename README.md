# pnp
So I keep seeing people using very heavy robust loss functions for pnp, instead of a well considered ransac loop. 
In my experience the latter wins in every case except very high inlier noise, something which doesnt actually happen in practice for pnp problems. 

This is a simple but fast pnp solver which for most problems does not require you to think at all. 
By default it will work well in most cases, but have a look in the parameters file for 
