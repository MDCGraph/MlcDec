//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_DATABUF_H
#define MLCDEC_DATABUF_H

#include "../Header.h"

#define PAGE_SIZE (1<<12)

// Page containing massive small/medium-sized data lists.
template<class T>
struct Page {
    T data[PAGE_SIZE]{};
    Page *next_page{nullptr};
};

// Buffer for small/medium-sized data.
template<class T>
struct MDataBuf {

    int cur_pos{0};
    Page<T> *head;
    Page<T> *cur_page;

    inline T *Allocate(size_t len = 1) {
        T *space;

        if (cur_pos + len > PAGE_SIZE) {
            auto page = new Page<T>();
            cur_page->next_page = page;
            cur_page = page;
            cur_pos = 0;
        }

        space = cur_page->data + cur_pos;
        cur_pos += len;
        return space;
    }

    MDataBuf() {
        head = new Page<T>();
        cur_page = head;
    }

    ~MDataBuf() = default;

    void Merge(MDataBuf<T> & d_buf) {
        cur_page->next_page = d_buf.head;
        cur_page = d_buf.cur_page;
        cur_pos = d_buf.cur_pos;
    }

    void Release() {
        cur_page = head;
        while (cur_page) {
            head = head->next_page;
            delete cur_page;
            cur_page = head;
        }
    }
};

// Page containing a single large-sized data list.
template<class T>
struct LPage {
    T *data{};
    LPage<T> *next_page{nullptr};

    LPage() = default;

    explicit LPage(size_t len) {
        data = new T[len];
    }

    ~LPage() {
        delete[] data;
    }
};

// Buffer for large-sized data.
template<class T>
struct LDataBuf {
    LPage<T> *head;
    LPage<T> *cur_page;

    inline T *Allocate(size_t len = 1) {
        auto page = new LPage<T>(len);
        cur_page->next_page = page;
        cur_page = page;
        return page->data;
    }
    
    LDataBuf() {
        head = new LPage<T>();
        cur_page = head;
    }

    ~LDataBuf() = default;

    void Merge(LDataBuf<T> & d_buf) {
        cur_page->next_page = d_buf.head;
        cur_page = d_buf.cur_page;
    }

    void Release() {
        cur_page = head;
        while (cur_page) {
            head = head->next_page;
            delete cur_page;
            cur_page = head;
        }
    }
};

template<class T>
struct DataBuf {

    MDataBuf<T> md_buf;
    LDataBuf<T> ld_buf;

    inline T *Allocate(size_t len) {
        return len >= PAGE_SIZE ? ld_buf.Allocate(len) : md_buf.Allocate(len);
    }

    DataBuf() = default;

    ~DataBuf() = default;

    void Merge(DataBuf<T> & d_buf) {
        md_buf.Merge(d_buf.md_buf);
        ld_buf.Merge(d_buf.ld_buf);
    }

    void Release() {
        md_buf.Release();
        ld_buf.Release();
    }
};


#endif //MLCDEC_DATABUF_H
