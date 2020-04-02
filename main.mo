model Main
   import Modelica.Constants.inf;
   import Pi=Modelica.Constants.pi;
   import pi=Modelica.Constants.pi;

   inner parameter Modelica.SIunits.Density rhoFluid = FP1.rhoFluid;
   inner parameter Modelica.SIunits.BulkModulus ElFluid = FP1.ElFluid;
   inner parameter Modelica.SIunits.KinematicViscosity nuFluid = FP1.nuFluid;
   public inner Maplesoft.Multibody.World world(gravityDir=Maplesoft.Multibody.Selectors.UnitVector.negY, gravityAcc=9.81);
   public Modelica.Blocks.Sources.Sine S1(amplitude=1.5, freqHz=.108333333333, phase=3.14159265359, offset=0, startTime=200) ;
   public Modelica.Mechanics.Translational.Sources.Position P1(useSupport=false, exact=true, f_crit=50) ;
   public Modelica.Mechanics.Translational.Components.Fixed F3(s0=0) ;
   public Maplesoft.Hydraulics.Basic.HydraulicCylinder HC1(A=0.113411494795e-2) ;
   public Maplesoft.Hydraulics.Basic.CheckValve CV2(Ropen=178.945936903, Gclosed=0.10e-11, Startclosed=false) ;
   public Maplesoft.Hydraulics.Basic.CheckValve CV1(Ropen=178.945936903, Gclosed=0.10e-11, Startclosed=true) ;
   public Maplesoft.Hydraulics.Basic.AtmosphericPressure AP1(P=.0) ;
   public Maplesoft.Hydraulics.Basic.CircularPipe CP1(D=.1003, L=1500.0, epsilon=0.15e-4, ReL=2000, ReT=4000, rho=rhoFluid, nu=nuFluid) ;
   public Maplesoft.Hydraulics.Basic.AtmosphericPressure AP2(P=100000.0) ;
   public inner Maplesoft.Hydraulics.FluidProperties FP1(rhoFluid=1000.0, ElFluid=800000000.0, nuFluid=0.2e-5) ;
   public MapleSimStandaloneSubsystem n0 ;
   public Maplesoft.Hydraulics.Basic.FluidInertia FI1(A=0.761764747263e-2, L=1500.0, rho=rhoFluid) ;
   public Maplesoft.Mechanical.Basic.TransFriction TF1(fs=702.75, fc=562.2, d=0, vs=.1, n=1, v0=0.1e-1) ;
   public MapleSimStandaloneSubsystem n1 ;
   public MapleSimStandaloneSubsystem n2 ;
   public MapleSimStandaloneSubsystem n3 ;
   public MapleSimStandaloneSubsystem n4 ;
   public MapleSimStandaloneSubsystem n5 ;
   public MapleSimStandaloneSubsystem n6 ;
   public MapleSimStandaloneSubsystem n7 ;
   public Maplesoft.Hydraulics.Basic.FixedPressure FP2(P=0.13234e8) ;
equation
   connect(S1.y, P1.s_ref) ;
   connect(CV2.portA, AP1.portA) ;
   connect(CV2.portB, CV1.portA) ;
   connect(F3.flange, HC1.flange_b) ;
   connect(HC1.portA, CV2.portB) ;
   connect(P1.flange, n0.b1) ;
   connect(CP1.portB, FI1.portA) ;
   connect(FI1.portB, AP2.portA) ;
   connect(HC1.flange_a, TF1.flange_b) ;
   connect(TF1.flange_a, F3.flange) ;
   connect(n1.b1, n0.a1) ;
   connect(n2.b1, n1.a1) ;
   connect(n3.b1, n2.a1) ;
   connect(n3.a1, n4.b1) ;
   connect(n4.a1, n5.b1) ;
   connect(n5.a1, n6.b1) ;
   connect(n6.a1, n7.b1) ;
   connect(n7.a1, HC1.flange_a) ;
   connect(FP2.portB, CP1.portA) ;
   connect(CV1.portB, FP2.portA) ;
end Main;
model VerticalPipe
    extends Maplesoft.Icons.CustomComponent;
    parameter Real z = 0 "z";
    parameter Real rho = 1000 "rho";
    parameter Real g = 9.81 "g";
    Real Pin;
    Real Qin;
    Real dP;
    Real Pout;
    Real Qout;
    Maplesoft.Hydraulics.Interfaces.Port_a portA;
    Maplesoft.Hydraulics.Interfaces.Port_a portB;
equation
    Qin = Qout;
    dP = rho * g * z;
    Pin - Pout = dP;
    portA.p = Pin;
    portA.q = Qin;
    portB.p = Pout;
    portB.q = -Qout;
end VerticalPipe;
model MapleSimStandaloneSubsystem
   import Modelica.Constants.inf;
   import Pi=Modelica.Constants.pi;
   import pi=Modelica.Constants.pi;

   public Modelica.Mechanics.Translational.Components.Mass M1(m=425.29310548, L=0) ;
   public Maplesoft.Sensors.SpeedSensor SS1(dimension="Velocity", fromUnit="m/s", toUnit="m/s", scale=1.0, shift=.0) ;
   public Modelica.Mechanics.Translational.Sources.Force F2(useSupport=false) ;
   public Modelica.Blocks.Math.BooleanToReal BTR2(realTrue=1.0, realFalse=.0) ;
   public Modelica.Blocks.Math.Gain G2(k=46.9936711616) ;
   public Modelica.Blocks.Logical.LessThreshold LT1(threshold=0) ;
   public Modelica.Mechanics.Translational.Sources.ConstantForce CF1(useSupport=false, f_constant=-3650.60969416) ;
   public Modelica.Blocks.Math.Product P3 ;
   public Modelica.Mechanics.Translational.Interfaces.Flange_b b1 ;
   public Modelica.Mechanics.Translational.Interfaces.Flange_a a1 ;
   public Maplesoft.Mechanical.Basic.TransFriction TF1(fs=0., fc=0., d=46.9936711616, vs=.1, n=1, v0=0.1e-1) ;
   public Modelica.Mechanics.Translational.Components.Fixed F1(s0=0) ;
   public Modelica.Blocks.Math.Gain G3(k=100) ;
   public Modelica.Blocks.Math.Tanh T1 ;
   public Modelica.Blocks.Math.Product P1 ;
   protected Modelica.Blocks.Interfaces.RealOutput RO_1 ;
   public Maplesoft.Sensors.ForceSensor FS1(dimension="Force", fromUnit="N", toUnit="N", scale=1.0, shift=.0) ;
   public Maplesoft.Sensors.ForceSensor FS2(dimension="Force", fromUnit="N", toUnit="N", scale=1.0, shift=.0) ;
   public Modelica.Blocks.Math.Gain G1(k=.0) ;
   public Modelica.Blocks.Math.Gain G4(k=.0) ;
   public Modelica.Blocks.Math.Add A1(k1=1, k2=1) ;
   public Modelica.Blocks.Math.Abs abs1_1 ;
   public Modelica.Blocks.Math.Add A2(k1=1, k2=-1) ;
   public Modelica.Mechanics.Translational.Components.SpringDamper SD1(s_nominal=0.10e-3, c=317552.185425, d=14849.8907517, s_rel0=0) ;
   //public outer Maplesoft.Hydraulics.FluidProperties FP1;
equation
   connect(M1.flange_a, SS1.F) ;
   connect(CF1.flange, M1.flange_a) ;
   connect(LT1.y, BTR2.u) ;
   connect(BTR2.y, G2.u) ;
   connect(F2.flange, M1.flange_a) ;
   connect(F1.flange, TF1.flange_a) ;
   connect(G3.y, T1.u) ;
   connect(P1.y, F2.f) ;
   connect(TF1.flange_b, M1.flange_a) ;
   connect(SS1.RO, LT1.u) ;
   connect(P3.u2, G2.y) ;
   connect(G3.u, RO_1) ;
   connect(SS1.RO, RO_1) ;
   connect(RO_1, P3.u1) ;
   connect(P1.u1, T1.y) ;
   connect(P3.y, A2.u1) ;
   connect(P1.u2, A2.y) ;
   connect(A2.u2, abs1_1.y) ;
   connect(A1.y, abs1_1.u) ;
   connect(G4.y, A1.u2) ;
   connect(G1.y, A1.u1) ;
   connect(FS1.RO, G1.u) ;
   connect(FS1.F1, b1) ;
   connect(FS2.F1, M1.flange_a) ;
   connect(FS2.RO, G4.u) ;
   connect(FS2.F2, a1) ;
   connect(SD1.flange_b, FS1.F2) ;
   connect(M1.flange_b, SD1.flange_a) ;
end MapleSimStandaloneSubsystem;
