//
// Created by janos on 9/14/20.
//

#pragma once

#include "Types.h"

#include <Magnum/Magnum.h>

#include <mutex>
#include <thread>
#include <condition_variable>
#include <queue>

extern "C" {

#include <libavcodec/avcodec.h>
#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>

}

namespace Phasefield {

namespace Mg = Magnum;

/**
 * aynchronous video saver
*/

class VideoSaver {
public:

    void startRecording(const char* path, Vector2i const& size);

    void endRecording();

    void appendFrame(Mg::Image2D&& frame);

    void encodeSingleFrame(Mg::Image2D const& frame);

    void work();

private:

    std::thread m_thread;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_finished = false;
    std::queue<Mg::Image2D> m_buffer;

    AVCodecContext* enc_ctx = nullptr;
    AVFrame* frame = nullptr;
    AVPacket* pkt = nullptr;
    AVCodecContext *c = nullptr;
    const AVCodec *codec;
    SwsContext *sws_ctx;
    FILE* outfile = nullptr;

    static constexpr uint8_t endcode[] = { 0, 0, 1, 0xb7 };

    size_t encodedFrames = 0;

    uint8_t *m_data[4];
    int m_line_size[4];
};

}


