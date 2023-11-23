//
// Created by ldd on 2023/3/6.
//

/*
 * The linear heap data structure is proposed by Chang et al. in:
 *
 * " Lijun Chang and Lu Qin. Cohesive Subgraph Computation over
 * Large Sparse Graphs. Springer Series in the Data Sciences, 2018 "
 *
 * Below is the array-based implementation, referencing the source code in
 * https://github.com/LijunChang/Cohesive_subgraph_book
 *
 */


#ifndef MLCDEC_INTLINEARHEAP_H
#define MLCDEC_INTLINEARHEAP_H


#include "../Header.h"

class IntLinearHeap {

public:
    IntLinearHeap() = default;

    ~IntLinearHeap() {
        if (release) {
            delete[] pre;
            delete[] next;
            delete[] heap;
        }

        if (release_bins) {
            delete[] bin;
        }
    }

    void Init(uint n_) {
        n = n_;

        pre = new uint[n];
        next = new uint[n];
        heap = new uint[n];
    }

    void Init(uint n_, uint *buf) {
        n = n_;

        pre = buf;
        next = buf + n;
        heap = buf + (n << 1);

        release = false;
    }


    void SetBin(uint max_bin_) {
        max_bin = max_bin_;
        bin = new uint[max_bin + 1];
        Clear();
    }

    void SetBin(uint max_bin_, uint *buf) {
        max_bin = max_bin_;
        bin = buf;
        Clear();

        release_bins = false;
    }

    void Copy(IntLinearHeap &ilh) {
        memcpy(pre, ilh.pre, n * sizeof(uint));
        memcpy(next, ilh.next, n * sizeof(uint));
        memcpy(heap, ilh.heap, n * sizeof(uint));

        memcpy(bin, ilh.bin, (max_bin + 1) * sizeof(uint));

        max_value = ilh.max_value;
        min_value = ilh.min_value;
    }

    inline void ClearBin(uint bid) {
        bin[bid] = -1;
    }

    inline void Insert(uint v, uint value) {
        uint head;

        head = bin[value];
        heap[v] = value;
        pre[v] = -1;
        next[v] = head;
        bin[value] = v;
        if (head != -1) pre[head] = v;

        if (value < min_value) min_value = value;
        if (value > max_value) max_value = value;
    }

    inline void Remove(uint v) {
        if (pre[v] == -1) { // v is the first element of the bin
            bin[heap[v]] = next[v];
            if (next[v] != -1) pre[next[v]] = -1;
        } else { // v is not the first element
            uint pv = pre[v];
            next[pv] = next[v];
            if (next[v] != -1) pre[next[v]] = pv;
        }
    }

    inline void Update(uint v, uint new_value) {
        Remove(v);
        Insert(v, new_value);
    }

    inline uint GetElements(uint value, uint *arr) {
        uint curr, i;

        i = 0;
        if (value >= min_value && value <= max_value) {
            curr = bin[value];
            while (curr != -1) {
                arr[i++] = curr;
                curr = next[curr];
            }
        }
        return i;
    }

    [[nodiscard]] uint GetMinValue() {
        Tighten();
        return min_value;
    }

    [[nodiscard]] uint GetMaxValue() {
        Tighten();
        return max_value;
    }

    uint GetValue(uint v) {
        return heap[v];
    }

    inline bool Empty() {
        Tighten();
        return min_value > max_value;
    }

    inline void Clear() {
        for (uint i = 0; i <= max_bin; i++) bin[i] = -1;
        min_value = max_bin;
        max_value = 0;
    }

private:
    inline void Tighten() {
        while (min_value <= max_value && bin[min_value] == -1) ++min_value;
        while (min_value <= max_value && bin[max_value] == -1) --max_value;
    }

    uint *pre{nullptr};
    uint *next{nullptr};
    uint *bin{nullptr};
    uint *heap{nullptr};

    uint n{};
    uint max_bin{};
    uint max_value{}, min_value{};

    bool release{true};
    bool release_bins{true};
};


#endif //MLCDEC_INTLINEARHEAP_H
