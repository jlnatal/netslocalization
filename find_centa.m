function centah=find_centa(state,positions,a)

theta_x=(positions(find(state>a),1))/5*2*pi;	%vector of x positions, mapped onto an angle
theta_y=(positions(find(state>a),2))/5*2*pi;	%vector of x positions, mapped onto an angle

squiggle_x=cos(theta_x)*5/(2*pi);	%weighted by the mass
curve_x=sin(theta_x)*5/(2*pi);		%weighted by the mass

squiggle_y=cos(theta_y)*5/(2*pi);		%weighted by the mass
curve_y=sin(theta_y)*5/(2*pi);		%weighted by the mass

squiggle_x_avg=1/(sum((find(state>a)),2))*sum((find(state>a)).*squiggle_x');
curve_x_avg=1/(sum((find(state>a)),2))*sum((find(state>a)).*curve_x');
squiggle_y_avg=1/(sum((find(state>a)),2))*sum((find(state>a)).*squiggle_y');
curve_y_avg=1/(sum((find(state>a)),2))*sum((find(state>a)).*curve_y');



theta_hat_x=atan2(-curve_x_avg,-squiggle_x_avg)+pi;
theta_hat_y=atan2(-curve_y_avg,-squiggle_y_avg)+pi;

x_com=5*theta_hat_x/(2*pi);
y_com=5*theta_hat_y/(2*pi);

centah=[x_com y_com];

