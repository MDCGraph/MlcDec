//
// Created by ldd on 2023/9/25.
//

#ifndef MLCDEC_CONSTRUCTS_H
#define MLCDEC_CONSTRUCTS_H


static uint log2_up(uint i) {
    uint a = 0;
    uint b = i - 1;
    while (b > 0) {
        b = b >> 1;
        a++;
    }
    return a;
}

static uint hash32(uint32_t a) {
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) ^ (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a + 0xd3a2646c) ^ (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) ^ (a >> 16);
    return a;
}


#endif //MLCDEC_CONSTRUCTS_H
