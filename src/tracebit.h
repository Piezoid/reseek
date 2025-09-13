#ifndef tracebit_h
#define tracebit_h

const byte TRACEBITS_DM = 0x01;
const byte TRACEBITS_IM = 0x02;
const byte TRACEBITS_MD = 0x04;
const byte TRACEBITS_MI = 0x08;
const byte TRACEBITS_SM = 0x10;
const byte TRACEBITS_DD = 0x20;  // Explicit D->D gap extension
const byte TRACEBITS_II = 0x40;  // Explicit I->I gap extension
const byte TRACEBITS_UNINIT = ~0x7f;

#endif // tracebit_h
