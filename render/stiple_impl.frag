//VTK::TCoord::Impl
int stip = 0x1 & (StipplePattern >> int(fract(stippleCoord)*16.0));
if (stip == 0)
   discard;