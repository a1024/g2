G2_R_R(logic_not){IDX; ASSIGN_R(xr[idx]==0);}
G2_R_C(logic_not){IDX; ASSIGN_R(!istrue_c(VEC2(x)));}
G2_R_Q(logic_not){IDX; ASSIGN_R(!istrue_q(VEC4(x)));}
DISC_R_O(logic_not){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_RR(logic_and){IDX; ASSIGN_R(xr[idx]&&yr[idx]);}
G2_R_RC(logic_and){IDX; ASSIGN_R(xr[idx]&&istrue_c(VEC2(y)));}
G2_R_RQ(logic_and){IDX; ASSIGN_R(xr[idx]&&istrue_q(VEC4(y)));}
G2_R_CR(logic_and){IDX; ASSIGN_R(istrue_c(VEC2(x))&&yr[idx]);}
G2_R_CC(logic_and){IDX; ASSIGN_R(istrue_c(VEC2(x))&&istrue_c(VEC2(y)));}
G2_R_CQ(logic_and){IDX; ASSIGN_R(istrue_c(VEC2(x))&&istrue_q(VEC4(y)));}
G2_R_QR(logic_and){IDX; ASSIGN_R(istrue_q(VEC4(x))&&yr[idx]);}
G2_R_QC(logic_and){IDX; ASSIGN_R(istrue_q(VEC4(x))&&istrue_c(VEC2(y)));}
G2_R_QQ(logic_and){IDX; ASSIGN_R(istrue_q(VEC4(x))&&istrue_q(VEC4(y)));}
DISC_R_O(logic_and){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_RR(logic_or){IDX; ASSIGN_R(xr[idx]||yr[idx]);}
G2_R_RC(logic_or){IDX; ASSIGN_R(xr[idx]||istrue_c(VEC2(y)));}
G2_R_RQ(logic_or){IDX; ASSIGN_R(xr[idx]||istrue_q(VEC4(y)));}
G2_R_CR(logic_or){IDX; ASSIGN_R(istrue_c(VEC2(x))||yr[idx]);}
G2_R_CC(logic_or){IDX; ASSIGN_R(istrue_c(VEC2(x))||istrue_c(VEC2(y)));}
G2_R_CQ(logic_or){IDX; ASSIGN_R(istrue_c(VEC2(x))||istrue_q(VEC4(y)));}
G2_R_QR(logic_or){IDX; ASSIGN_R(istrue_q(VEC4(x))||yr[idx]);}
G2_R_QC(logic_or){IDX; ASSIGN_R(istrue_q(VEC4(x))||istrue_c(VEC2(y)));}
G2_R_QQ(logic_or){IDX; ASSIGN_R(istrue_q(VEC4(x))||istrue_q(VEC4(y)));}
DISC_R_O(logic_or){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_RR(logic_xor){IDX; ASSIGN_R((xr[idx]!=0)^(yr[idx]!=0));}
G2_R_RC(logic_xor){IDX; ASSIGN_R((xr[idx]!=0)^istrue_c(VEC2(y)));}
G2_R_RQ(logic_xor){IDX; ASSIGN_R((xr[idx]!=0)^istrue_q(VEC4(y)));}
G2_R_CR(logic_xor){IDX; ASSIGN_R(istrue_c(VEC2(x))^(yr[idx]!=0));}
G2_R_CC(logic_xor){IDX; ASSIGN_R(istrue_c(VEC2(x))^istrue_c(VEC2(y)));}
G2_R_CQ(logic_xor){IDX; ASSIGN_R(istrue_c(VEC2(x))^istrue_q(VEC4(y)));}
G2_R_QR(logic_xor){IDX; ASSIGN_R(istrue_q(VEC4(x))^(yr[idx]!=0));}
G2_R_QC(logic_xor){IDX; ASSIGN_R(istrue_q(VEC4(x))^istrue_c(VEC2(y)));}
G2_R_QQ(logic_xor){IDX; ASSIGN_R(istrue_q(VEC4(x))^istrue_q(VEC4(y)));}
DISC_R_O(logic_xor){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}