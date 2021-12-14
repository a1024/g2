G2_R_R(logic_equal){IDX; ASSIGN_R(xr[idx]==0);}
G2_R_C(logic_equal){IDX; ASSIGN_R(!istrue_c(VEC2(x)));}
G2_R_Q(logic_equal){IDX; ASSIGN_R(!istrue_q(VEC4(x)));}
G2_R_RR(logic_equal){IDX; ASSIGN_R(xr[idx]==yr[idx]);}
G2_R_RC(logic_equal){IDX; ASSIGN_R(equal_rc(xr[idx], VEC2(y)));}
G2_R_RQ(logic_equal){IDX; ASSIGN_R(equal_rq(xr[idx], VEC4(y)));}
G2_R_CR(logic_equal){IDX; ASSIGN_R(equal_cr(VEC2(x), yr[idx]));}
G2_R_CC(logic_equal){IDX; ASSIGN_R(equal_cc(VEC2(x), VEC2(y)));}
G2_R_CQ(logic_equal){IDX; ASSIGN_R(equal_cq(VEC2(x), VEC4(y)));}
G2_R_QR(logic_equal){IDX; ASSIGN_R(equal_qr(VEC4(x), yr[idx]));}
G2_R_QC(logic_equal){IDX; ASSIGN_R(equal_qc(VEC4(x), VEC2(y)));}
G2_R_QQ(logic_equal){IDX; ASSIGN_R(equal_qq(VEC4(x), VEC4(y)));}
DISC_R_O(logic_equal){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_R(logic_not_equal){IDX; ASSIGN_R(xr[idx]!=0);}
G2_R_C(logic_not_equal){IDX; ASSIGN_R(istrue_c(VEC2(x)));}
G2_R_Q(logic_not_equal){IDX; ASSIGN_R(istrue_q(VEC4(x)));}
G2_R_RR(logic_not_equal){IDX; ASSIGN_R(xr[idx]!=yr[idx]);}
G2_R_RC(logic_not_equal){IDX; ASSIGN_R(!equal_rc(xr[idx], VEC2(y)));}
G2_R_RQ(logic_not_equal){IDX; ASSIGN_R(!equal_rq(xr[idx], VEC4(y)));}
G2_R_CR(logic_not_equal){IDX; ASSIGN_R(!equal_cr(VEC2(x), yr[idx]));}
G2_R_CC(logic_not_equal){IDX; ASSIGN_R(!equal_cc(VEC2(x), VEC2(y)));}
G2_R_CQ(logic_not_equal){IDX; ASSIGN_R(!equal_cq(VEC2(x), VEC4(y)));}
G2_R_QR(logic_not_equal){IDX; ASSIGN_R(!equal_qr(VEC4(x), yr[idx]));}
G2_R_QC(logic_not_equal){IDX; ASSIGN_R(!equal_qc(VEC4(x), VEC2(y)));}
G2_R_QQ(logic_not_equal){IDX; ASSIGN_R(!equal_qq(VEC4(x), VEC4(y)));}
DISC_R_O(logic_not_equal){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_R(logic_less_l){IDX; ASSIGN_R(0<xr[idx]);}
G2_R_C(logic_less_l){IDX; ASSIGN_R(0<xr[idx]);}
G2_R_Q(logic_less_l){IDX; ASSIGN_R(0<xr[idx]);}
G2_R_R(logic_less_r){IDX; ASSIGN_R(xr[idx]<0);}
G2_R_C(logic_less_r){IDX; ASSIGN_R(xr[idx]<0);}
G2_R_Q(logic_less_r){IDX; ASSIGN_R(xr[idx]<0);}
G2_R_RR(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_RC(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_RQ(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_CR(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_CC(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_CQ(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_QR(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_QC(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
G2_R_QQ(logic_less){IDX; ASSIGN_R(xr[idx]<yr[idx]);}
DISC_R_O(logic_less){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}