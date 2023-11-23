//
// Created by ldd on 2023/10/28.
//

#ifndef MLCDEC_MLCTREEBUILDOPTION_H
#define MLCDEC_MLCTREEBUILDOPTION_H

enum RUN_OPTION {
    IP, IP_MERGE, IP_START, P_OPT, SERIAL
};

static string RunOpt2Str(RUN_OPTION opt) {
    if (opt == IP) return "ip";
    else if (opt == IP_MERGE) return "ipm";
    else if (opt == IP_START) return "ips";
    else if (opt == P_OPT) return "opt";
    else if (opt == SERIAL) return "serial";
    else {
        cerr << "Invalid running option." << endl;
        exit(-1);
    }
}


#endif //MLCDEC_MLCTREEBUILDOPTION_H
