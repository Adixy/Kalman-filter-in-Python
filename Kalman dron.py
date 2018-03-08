import math
class Kalman:
    def __init__(self,init_cond=[0,0,0,0,0,0],m=3.7,Ix=0.04,Iy=0.04,Iz=0.05):
        self.init_cond=init_cond
        self.m=m
        self.ix=Ix
        self.iy=Iy
        self.iz=Iz
        self.fi_p1=0
        self.fi_p2=0
        self.fi_p3=0
        self.fi_p4=0
        self.fi_p5=0
        self.fi_p6=0
        self.fi_p7=0
        self.fi_p8=0
        self.fi_p9=0
        self.th_p1=0
        self.th_p2=0
        self.th_p3=0
        self.th_p4=0
        self.th_p5=0
        self.th_p6=0
        self.th_p7=0
        self.th_p8=0
        self.th_p9=0
        self.psi_p1=0
        self.psi_p2=0
        self.psi_p3=0
        self.psi_p4=0
        self.psi_p5=0
        self.psi_p6=0
        self.psi_p7=0
        self.psi_p8=0
        self.psi_p9=0
    def update_rot(self,acc,gyro,mag=[0,0,0],dt,tau):
        d2t=dt*dt
        d3t=dt*dt*dt
        d4t=dt*dt*dt*dt
        tau_x=tau[0]
        tau_y=tau[1]
        tau_z=tau[2]
        
        self.d2fi=(tau_x/self.ix)+self.d2fi
        self.dfi=dt*self.d2fi+self.dfi
        self.fi=(d2t/2)*self.d2fi+dt*self.dfi+self.fi
        
        self.fi_p1=self_fi_q1+self.fi_p1
        self.fi_p2=self.fi_q2+dt*self.fi_p1+self.fi_p2
        self.fi_p3=self.fi_q3+(d2t/2)*self.fi_p1+dt*self.fi_p2+self.fi_p3
        self.fi_p4=self.fi_q4+self.fi_p1*dt+self.fi_p4
        self.fi_p5=self.fi_q5+d2t*self.fi_p1+self.fi_p4+dt*self.fi_p2+self.fi_p5
        self.fi_p6=self.fi_q6+(d3t/2)*self.fi_p1+(d2t/2)*self.fi_p4+d2t*self.fi_p2+dt*self.fi_p5+dt*self.fi_p3+self.fi_p6
        self.fi_p7=self.fi_q7+(d2t/2)*self.fi_p1+dt*self.fi_p4+self.fi_p7
        self.fi_p8=self.fi_q8+(d3t/2)*self.fi_p1+d2t*self.fi_p4+dt*self.fi_p7+(d2t/2)*self.fi_p2+dt*self.fi_p5+self.fi_p8
        self.fi_p9=self.fi_q9+(d4t/4)*self.fi_p1+(d3t/2)*self.fi_p4+(d2t/2)*self.fi_p7+(d3t/2)*self.fi_p2+d2t*self.fi_p5+dt*self.fi_p8+(d2t/2)*self.fi_p3+dt*self.fi_p6+self.fi_p9

        T1=self.fi_p1+self.fi_r1
        T2=self.fi_p2+self.fi_r2
        T3=self.fi_p3+self.fi_r3
        T4=self.fi_p4+self.fi_r4
        T5=self.fi_p5+self.fi_r5
        T6=self.fi_p6+self.fi_r6
        T7=self.fi_p7+self.fi_r7
        T8=self.fi_p8+self.fi_r8
        T9=self.fi_p9+self.fi_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.fi_p1*i1+self.fi_p2*i4+self.fi_p3*i7
            K2=self.fi_p1*i2+self.fi_p2*i5+self.fi_p3*i8
            K3=self.fi_p1*i3+self.fi_p2*i6+self.fi_p3*i9
            K4=self.fi_p4*i1+self.fi_p5*i4+self.fi_p6*i7
            K5=self.fi_p4*i2+self.fi_p5*i5+self.fi_p6*i8
            K6=self.fi_p4*i3+self.fi_p5*i6+self.fi_p6*i9
            K7=self.fi_p7*i1+self.fi_p8*i4+self.fi_p9*i7
            K8=self.fi_p7*i2+self.fi_p8*i5+self.fi_p9*i8
            K9=self.fi_p7*i3+self.fi_p8*i6+self.fi_p9*i9
            
        else:
            print('Error 1.1, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2fi=self.d2fi+K1*(zfi1-self.d2fi)+K2*(zfi2-self.d1fi)+K3*(zfi3-self.fi)
        self.d1fi=self.d1fi+K4*(zfi1-self.d2fi)+K5*(zfi2-self.d1fi)+K6*(zfi3-self.fi)
        self.fi=self.fi+K7*(zfi1-self.d2fi)+K8*(zfi2-self.d1fi)+K9*(zfi3-self.fi)

        self.fi_p1=self.fi_p1*(1-K1)-self.fi_p4*K2-self.fi_p7*K3
        self.fi_p2=self.fi_p2*(1-K1)-self.fi_p5*K2-self.fi_p8*K3
        self.fi_p3=self.fi_p3*(1-K1)-self.fi_p6*K2-self.fi_p9*K3
        self.fi_p4=-self.fi_p1*K4+self.fi_p4*(1-K5)-self.fi_p7*K6
        self.fi_p5=-self.fi_p2*K4+self.fi_p5*(1-K5)-self.fi_p8*K6
        self.fi_p6=-self.fi_p3*K4+self.fi_p6*(1-K5)-self.fi_p9*K6
        self.fi_p7=-self.fi_p1*K7-self.fi_p4*K8+self.fi_p7*(1-K9)
        self.fi_p8=-self.fi_p2*K7-self.fi_p5*K8+self.fi_p8*(1-K9)
        self.fi_p9=-self.fi_p3*K7-self.fi_p6*K8+self.fi_p9*(1-K9)
#######################################################################
        self.d2th=(tau_y/self.iy)+self.d2th
        self.dth=dt*self.d2th+self.dth
        self.th=(d2t/2)*self.d2th+dt*self.dth+self.th
        
        self.th_p1=self_th_q1+self.th_p1
        self.th_p2=self.th_q2+dt*self.th_p1+self.th_p2
        self.th_p3=self.th_q3+(d2t/2)*self.th_p1+dt*self.th_p2+self.th_p3
        self.th_p4=self.th_q4+self.th_p1*dt+self.th_p4
        self.th_p5=self.th_q5+d2t*self.th_p1+self.th_p4+dt*self.th_p2+self.th_p5
        self.th_p6=self.th_q6+(d3t/2)*self.th_p1+(d2t/2)*self.th_p4+d2t*self.th_p2+dt*self.th_p5+dt*self.th_p3+self.th_p6
        self.th_p7=self.th_q7+(d2t/2)*self.th_p1+dt*self.th_p4+self.th_p7
        self.th_p8=self.th_q8+(d3t/2)*self.th_p1+d2t*self.th_p4+dt*self.th_p7+(d2t/2)*self.th_p2+dt*self.th_p5+self.th_p8
        self.th_p9=self.th_q9+(d4t/4)*self.th_p1+(d3t/2)*self.th_p4+(d2t/2)*self.th_p7+(d3t/2)*self.th_p2+d2t*self.th_p5+dt*self.th_p8+(d2t/2)*self.th_p3+dt*self.th_p6+self.th_p9

        T1=self.th_p1+self.th_r1
        T2=self.th_p2+self.th_r2
        T3=self.th_p3+self.th_r3
        T4=self.th_p4+self.th_r4
        T5=self.th_p5+self.th_r5
        T6=self.th_p6+self.th_r6
        T7=self.th_p7+self.th_r7
        T8=self.th_p8+self.th_r8
        T9=self.th_p9+self.th_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.th_p1*i1+self.th_p2*i4+self.th_p3*i7
            K2=self.th_p1*i2+self.th_p2*i5+self.th_p3*i8
            K3=self.th_p1*i3+self.th_p2*i6+self.th_p3*i9
            K4=self.th_p4*i1+self.th_p5*i4+self.th_p6*i7
            K5=self.th_p4*i2+self.th_p5*i5+self.th_p6*i8
            K6=self.th_p4*i3+self.th_p5*i6+self.th_p6*i9
            K7=self.th_p7*i1+self.th_p8*i4+self.th_p9*i7
            K8=self.th_p7*i2+self.th_p8*i5+self.th_p9*i8
            K9=self.th_p7*i3+self.th_p8*i6+self.th_p9*i9
            
        else:
            print('Error 1.2, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2th=self.d2th+K1*(zth1-self.d2th)+K2*(zth2-self.d1th)+K3*(zth3-self.th)
        self.d1th=self.d1th+K4*(zth1-self.d2th)+K5*(zth2-self.d1th)+K6*(zth3-self.th)
        self.th=self.th+K7*(zth1-self.d2th)+K8*(zth2-self.d1th)+K9*(zth3-self.th)

        self.th_p1=self.th_p1*(1-K1)-self.th_p4*K2-self.th_p7*K3
        self.th_p2=self.th_p2*(1-K1)-self.th_p5*K2-self.th_p8*K3
        self.th_p3=self.th_p3*(1-K1)-self.th_p6*K2-self.th_p9*K3
        self.th_p4=-self.th_p1*K4+self.th_p4*(1-K5)-self.th_p7*K6
        self.th_p5=-self.th_p2*K4+self.th_p5*(1-K5)-self.th_p8*K6
        self.th_p6=-self.th_p3*K4+self.th_p6*(1-K5)-self.th_p9*K6
        self.th_p7=-self.th_p1*K7-self.th_p4*K8+self.th_p7*(1-K9)
        self.th_p8=-self.th_p2*K7-self.th_p5*K8+self.th_p8*(1-K9)
        self.th_p9=-self.th_p3*K7-self.th_p6*K8+self.th_p9*(1-K9)
##########################################################################################
        self.d2psi=(tau_z/self.iz)+self.d2psi
        self.dpsi=dt*self.d2psi+self.dpsi
        self.psi=(d2t/2)*self.d2psi+dt*self.dpsi+self.psi
        
        self.psi_p1=self_psi_q1+self.psi_p1
        self.psi_p2=self.psi_q2+dt*self.psi_p1+self.psi_p2
        self.psi_p3=self.psi_q3+(d2t/2)*self.psi_p1+dt*self.psi_p2+self.psi_p3
        self.psi_p4=self.psi_q4+self.psi_p1*dt+self.psi_p4
        self.psi_p5=self.psi_q5+d2t*self.psi_p1+self.psi_p4+dt*self.psi_p2+self.psi_p5
        self.psi_p6=self.psi_q6+(d3t/2)*self.psi_p1+(d2t/2)*self.psi_p4+d2t*self.psi_p2+dt*self.psi_p5+dt*self.psi_p3+self.psi_p6
        self.psi_p7=self.psi_q7+(d2t/2)*self.psi_p1+dt*self.psi_p4+self.psi_p7
        self.psi_p8=self.psi_q8+(d3t/2)*self.psi_p1+d2t*self.psi_p4+dt*self.psi_p7+(d2t/2)*self.psi_p2+dt*self.psi_p5+self.psi_p8
        self.psi_p9=self.psi_q9+(d4t/4)*self.psi_p1+(d3t/2)*self.psi_p4+(d2t/2)*self.psi_p7+(d3t/2)*self.psi_p2+d2t*self.psi_p5+dt*self.psi_p8+(d2t/2)*self.psi_p3+dt*self.psi_p6+self.psi_p9

        T1=self.psi_p1+self.psi_r1
        T2=self.psi_p2+self.psi_r2
        T3=self.psi_p3+self.psi_r3
        T4=self.psi_p4+self.psi_r4
        T5=self.psi_p5+self.psi_r5
        T6=self.psi_p6+self.psi_r6
        T7=self.psi_p7+self.psi_r7
        T8=self.psi_p8+self.psi_r8
        T9=self.psi_p9+self.psi_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.psi_p1*i1+self.psi_p2*i4+self.psi_p3*i7
            K2=self.psi_p1*i2+self.psi_p2*i5+self.psi_p3*i8
            K3=self.psi_p1*i3+self.psi_p2*i6+self.psi_p3*i9
            K4=self.psi_p4*i1+self.psi_p5*i4+self.psi_p6*i7
            K5=self.psi_p4*i2+self.psi_p5*i5+self.psi_p6*i8
            K6=self.psi_p4*i3+self.psi_p5*i6+self.psi_p6*i9
            K7=self.psi_p7*i1+self.psi_p8*i4+self.psi_p9*i7
            K8=self.psi_p7*i2+self.psi_p8*i5+self.psi_p9*i8
            K9=self.psi_p7*i3+self.psi_p8*i6+self.psi_p9*i9
            
        else:
            print('Error 1.3, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2psi=self.d2psi+K1*(zpsi1-self.d2psi)+K2*(zpsi2-self.d1psi)+K3*(zpsi3-self.psi)
        self.d1psi=self.d1psi+K4*(zpsi1-self.d2psi)+K5*(zpsi2-self.d1psi)+K6*(zpsi3-self.psi)
        self.psi=self.psi+K7*(zpsi1-self.d2psi)+K8*(zpsi2-self.d1psi)+K9*(zpsi3-self.psi)

        self.psi_p1=self.psi_p1*(1-K1)-self.psi_p4*K2-self.psi_p7*K3
        self.psi_p2=self.psi_p2*(1-K1)-self.psi_p5*K2-self.psi_p8*K3
        self.psi_p3=self.psi_p3*(1-K1)-self.psi_p6*K2-self.psi_p9*K3
        self.psi_p4=-self.psi_p1*K4+self.psi_p4*(1-K5)-self.psi_p7*K6
        self.psi_p5=-self.psi_p2*K4+self.psi_p5*(1-K5)-self.psi_p8*K6
        self.psi_p6=-self.psi_p3*K4+self.psi_p6*(1-K5)-self.psi_p9*K6
        self.psi_p7=-self.psi_p1*K7-self.psi_p4*K8+self.psi_p7*(1-K9)
        self.psi_p8=-self.psi_p2*K7-self.psi_p5*K8+self.psi_p8*(1-K9)
        self.psi_p9=-self.psi_p3*K7-self.psi_p6*K8+self.psi_p9*(1-K9)
    def update_trans(self,acc,gyro,mag=[0,0,0],dt,Fm):
        d2t=dt*dt
        d3t=dt*dt*dt
        d4t=dt*dt*dt*dt
        
        self.d2x=(Fm/self.m)*math.sin(math.radians(self.th))+self.d2x
        self.dx=dt*self.d2x+self.dx
        self.x=(d2t/2)*self.d2x+dt*self.dx+self.x
        
        self.x_p1=self_x_q1+self.x_p1
        self.x_p2=self.x_q2+dt*self.x_p1+self.x_p2
        self.x_p3=self.x_q3+(d2t/2)*self.x_p1+dt*self.x_p2+self.x_p3
        self.x_p4=self.x_q4+self.x_p1*dt+self.x_p4
        self.x_p5=self.x_q5+d2t*self.x_p1+self.x_p4+dt*self.x_p2+self.x_p5
        self.x_p6=self.x_q6+(d3t/2)*self.x_p1+(d2t/2)*self.x_p4+d2t*self.x_p2+dt*self.x_p5+dt*self.x_p3+self.x_p6
        self.x_p7=self.x_q7+(d2t/2)*self.x_p1+dt*self.x_p4+self.x_p7
        self.x_p8=self.x_q8+(d3t/2)*self.x_p1+d2t*self.x_p4+dt*self.x_p7+(d2t/2)*self.x_p2+dt*self.x_p5+self.x_p8
        self.x_p9=self.x_q9+(d4t/4)*self.x_p1+(d3t/2)*self.x_p4+(d2t/2)*self.x_p7+(d3t/2)*self.x_p2+d2t*self.x_p5+dt*self.x_p8+(d2t/2)*self.x_p3+dt*self.x_p6+self.x_p9

        T1=self.x_p1+self.x_r1
        T2=self.x_p2+self.x_r2
        T3=self.x_p3+self.x_r3
        T4=self.x_p4+self.x_r4
        T5=self.x_p5+self.x_r5
        T6=self.x_p6+self.x_r6
        T7=self.x_p7+self.x_r7
        T8=self.x_p8+self.x_r8
        T9=self.x_p9+self.x_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.x_p1*i1+self.x_p2*i4+self.x_p3*i7
            K2=self.x_p1*i2+self.x_p2*i5+self.x_p3*i8
            K3=self.x_p1*i3+self.x_p2*i6+self.x_p3*i9
            K4=self.x_p4*i1+self.x_p5*i4+self.x_p6*i7
            K5=self.x_p4*i2+self.x_p5*i5+self.x_p6*i8
            K6=self.x_p4*i3+self.x_p5*i6+self.x_p6*i9
            K7=self.x_p7*i1+self.x_p8*i4+self.x_p9*i7
            K8=self.x_p7*i2+self.x_p8*i5+self.x_p9*i8
            K9=self.x_p7*i3+self.x_p8*i6+self.x_p9*i9
            
        else:
            print('Error 2.1, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2x=self.d2x+K1*(zx1-self.d2x)+K2*(zx2-self.d1x)+K3*(zx3-self.x)
        self.d1x=self.d1x+K4*(zx1-self.d2x)+K5*(zx2-self.d1x)+K6*(zx3-self.x)
        self.x=self.x+K7*(zx1-self.d2x)+K8*(zx2-self.d1x)+K9*(zx3-self.x)

        self.x_p1=self.x_p1*(1-K1)-self.x_p4*K2-self.x_p7*K3
        self.x_p2=self.x_p2*(1-K1)-self.x_p5*K2-self.x_p8*K3
        self.x_p3=self.x_p3*(1-K1)-self.x_p6*K2-self.x_p9*K3
        self.x_p4=-self.x_p1*K4+self.x_p4*(1-K5)-self.x_p7*K6
        self.x_p5=-self.x_p2*K4+self.x_p5*(1-K5)-self.x_p8*K6
        self.x_p6=-self.x_p3*K4+self.x_p6*(1-K5)-self.x_p9*K6
        self.x_p7=-self.x_p1*K7-self.x_p4*K8+self.x_p7*(1-K9)
        self.x_p8=-self.x_p2*K7-self.x_p5*K8+self.x_p8*(1-K9)
        self.x_p9=-self.x_p3*K7-self.x_p6*K8+self.x_p9*(1-K9)
#######################################################################
        self.d2y=-(Fm/self.m)*math.sin(math.radians(self.fi))+self.d2y
        self.dy=dt*self.d2y+self.dy
        self.y=(d2t/2)*self.d2y+dt*self.dy+self.y
        
        self.y_p1=self_y_q1+self.y_p1
        self.y_p2=self.y_q2+dt*self.y_p1+self.y_p2
        self.y_p3=self.y_q3+(d2t/2)*self.y_p1+dt*self.y_p2+self.y_p3
        self.y_p4=self.y_q4+self.y_p1*dt+self.y_p4
        self.y_p5=self.y_q5+d2t*self.y_p1+self.y_p4+dt*self.y_p2+self.y_p5
        self.y_p6=self.y_q6+(d3t/2)*self.y_p1+(d2t/2)*self.y_p4+d2t*self.y_p2+dt*self.y_p5+dt*self.y_p3+self.y_p6
        self.y_p7=self.y_q7+(d2t/2)*self.y_p1+dt*self.y_p4+self.y_p7
        self.y_p8=self.y_q8+(d3t/2)*self.y_p1+d2t*self.y_p4+dt*self.y_p7+(d2t/2)*self.y_p2+dt*self.y_p5+self.y_p8
        self.y_p9=self.y_q9+(d4t/4)*self.y_p1+(d3t/2)*self.y_p4+(d2t/2)*self.y_p7+(d3t/2)*self.y_p2+d2t*self.y_p5+dt*self.y_p8+(d2t/2)*self.y_p3+dt*self.y_p6+self.y_p9

        T1=self.y_p1+self.y_r1
        T2=self.y_p2+self.y_r2
        T3=self.y_p3+self.y_r3
        T4=self.y_p4+self.y_r4
        T5=self.y_p5+self.y_r5
        T6=self.y_p6+self.y_r6
        T7=self.y_p7+self.y_r7
        T8=self.y_p8+self.y_r8
        T9=self.y_p9+self.y_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.y_p1*i1+self.y_p2*i4+self.y_p3*i7
            K2=self.y_p1*i2+self.y_p2*i5+self.y_p3*i8
            K3=self.y_p1*i3+self.y_p2*i6+self.y_p3*i9
            K4=self.y_p4*i1+self.y_p5*i4+self.y_p6*i7
            K5=self.y_p4*i2+self.y_p5*i5+self.y_p6*i8
            K6=self.y_p4*i3+self.y_p5*i6+self.y_p6*i9
            K7=self.y_p7*i1+self.y_p8*i4+self.y_p9*i7
            K8=self.y_p7*i2+self.y_p8*i5+self.y_p9*i8
            K9=self.y_p7*i3+self.y_p8*i6+self.y_p9*i9
            
        else:
            print('Error 2.2, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2y=self.d2y+K1*(zy1-self.d2y)+K2*(zy2-self.d1y)+K3*(zy3-self.y)
        self.d1y=self.d1y+K4*(zy1-self.d2y)+K5*(zy2-self.d1y)+K6*(zy3-self.y)
        self.y=self.y+K7*(zy1-self.d2y)+K8*(zy2-self.d1y)+K9*(zy3-self.y)

        self.y_p1=self.y_p1*(1-K1)-self.y_p4*K2-self.y_p7*K3
        self.y_p2=self.y_p2*(1-K1)-self.y_p5*K2-self.y_p8*K3
        self.y_p3=self.y_p3*(1-K1)-self.y_p6*K2-self.y_p9*K3
        self.y_p4=-self.y_p1*K4+self.y_p4*(1-K5)-self.y_p7*K6
        self.y_p5=-self.y_p2*K4+self.y_p5*(1-K5)-self.y_p8*K6
        self.y_p6=-self.y_p3*K4+self.y_p6*(1-K5)-self.y_p9*K6
        self.y_p7=-self.y_p1*K7-self.y_p4*K8+self.y_p7*(1-K9)
        self.y_p8=-self.y_p2*K7-self.y_p5*K8+self.y_p8*(1-K9)
        self.y_p9=-self.y_p3*K7-self.y_p6*K8+self.y_p9*(1-K9)
##########################################################################################
        self.d2z=(Fm/self.m)*math.cos(math.radians(self.fi))*math.cos(math.radians(self.th))+self.d2z
        self.dz=dt*self.d2z+self.dz
        self.z=(d2t/2)*self.d2z+dt*self.dz+self.z
        
        self.z_p1=self_z_q1+self.z_p1
        self.z_p2=self.z_q2+dt*self.z_p1+self.z_p2
        self.z_p3=self.z_q3+(d2t/2)*self.z_p1+dt*self.z_p2+self.z_p3
        self.z_p4=self.z_q4+self.z_p1*dt+self.z_p4
        self.z_p5=self.z_q5+d2t*self.z_p1+self.z_p4+dt*self.z_p2+self.z_p5
        self.z_p6=self.z_q6+(d3t/2)*self.z_p1+(d2t/2)*self.z_p4+d2t*self.z_p2+dt*self.z_p5+dt*self.z_p3+self.z_p6
        self.z_p7=self.z_q7+(d2t/2)*self.z_p1+dt*self.z_p4+self.z_p7
        self.z_p8=self.z_q8+(d3t/2)*self.z_p1+d2t*self.z_p4+dt*self.z_p7+(d2t/2)*self.z_p2+dt*self.z_p5+self.z_p8
        self.z_p9=self.z_q9+(d4t/4)*self.z_p1+(d3t/2)*self.z_p4+(d2t/2)*self.z_p7+(d3t/2)*self.z_p2+d2t*self.z_p5+dt*self.z_p8+(d2t/2)*self.z_p3+dt*self.z_p6+self.z_p9

        T1=self.z_p1+self.z_r1
        T2=self.z_p2+self.z_r2
        T3=self.z_p3+self.z_r3
        T4=self.z_p4+self.z_r4
        T5=self.z_p5+self.z_r5
        T6=self.z_p6+self.z_r6
        T7=self.z_p7+self.z_r7
        T8=self.z_p8+self.z_r8
        T9=self.z_p9+self.z_r9

        D1=T1*(T5*T9-T6*T8)
        D1=T2*(T4*T9-T6*T7)
        D1=T3*(T4*T8-T5*T7)
        D=D1+D3-D2
        if D!=0:
            i1=(T5*T9-T6*T8)/D
            i2=(T6*T7-T4*T9)/D
            i3=(T4*T8-T5*T7)/D
            i4=(T3*T8-T2*T9)/D
            i5=(T1*T9-T3*T7)/D
            i6=(T2*T7-T1*T8)/D
            i7=(T2*T6-T3*T5)/D
            i8=(T3*T4-T1*T6)/D
            i9=(T1*T5-T2*T4)/D

            K1=self.z_p1*i1+self.z_p2*i4+self.z_p3*i7
            K2=self.z_p1*i2+self.z_p2*i5+self.z_p3*i8
            K3=self.z_p1*i3+self.z_p2*i6+self.z_p3*i9
            K4=self.z_p4*i1+self.z_p5*i4+self.z_p6*i7
            K5=self.z_p4*i2+self.z_p5*i5+self.z_p6*i8
            K6=self.z_p4*i3+self.z_p5*i6+self.z_p6*i9
            K7=self.z_p7*i1+self.z_p8*i4+self.z_p9*i7
            K8=self.z_p7*i2+self.z_p8*i5+self.z_p9*i8
            K9=self.z_p7*i3+self.z_p8*i6+self.z_p9*i9
            
        else:
            print('Error 2.3, determinant is zero!')
            K1=1
            K2=1
            K3=1
            K4=1
            K5=1
            K6=1
            K7=1
            K8=1
            K9=1
        self.d2z=self.d2z+K1*(zz1-self.d2z)+K2*(zz2-self.d1z)+K3*(zz3-self.z)
        self.d1z=self.d1z+K4*(zz1-self.d2z)+K5*(zz2-self.d1z)+K6*(zz3-self.z)
        self.z=self.z+K7*(zz1-self.d2z)+K8*(zz2-self.d1z)+K9*(zz3-self.z)

        self.z_p1=self.z_p1*(1-K1)-self.z_p4*K2-self.z_p7*K3
        self.z_p2=self.z_p2*(1-K1)-self.z_p5*K2-self.z_p8*K3
        self.z_p3=self.z_p3*(1-K1)-self.z_p6*K2-self.z_p9*K3
        self.z_p4=-self.z_p1*K4+self.z_p4*(1-K5)-self.z_p7*K6
        self.z_p5=-self.z_p2*K4+self.z_p5*(1-K5)-self.z_p8*K6
        self.z_p6=-self.z_p3*K4+self.z_p6*(1-K5)-self.z_p9*K6
        self.z_p7=-self.z_p1*K7-self.z_p4*K8+self.z_p7*(1-K9)
        self.z_p8=-self.z_p2*K7-self.z_p5*K8+self.z_p8*(1-K9)
        self.z_p9=-self.z_p3*K7-self.z_p6*K8+self.z_p9*(1-K9)
