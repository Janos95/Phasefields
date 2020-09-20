//
// Created by janos on 9/14/20.
//

#include "VideoSaver.h"
#include <Magnum/Image.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/Math/Color.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Corrade/Containers/Optional.h>

namespace Phasefield {

namespace {

void encode(AVCodecContext *enc_ctx, AVFrame *frame, AVPacket *pkt, FILE *outfile) {
    int ret;

    /* send the frame to the encoder */

    ret = avcodec_send_frame(enc_ctx, frame);
    if (ret < 0) {
        fprintf(stderr, "Error sending a frame for encoding\n");
        exit(1);
    }

    while (ret >= 0) {
        ret = avcodec_receive_packet(enc_ctx, pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            return;
        else if (ret < 0) {
            fprintf(stderr, "Error during encoding\n");
            exit(1);
        }

        fwrite(pkt->data, 1, pkt->size, outfile);
        av_packet_unref(pkt);
    }
}

}

void VideoSaver::endRecording() {
    {
        std::unique_lock lock(m_mutex);
        m_finished = true;
        m_cv.notify_all();
    }

    m_thread.join();

    /* flush the encoder */
    encode(c, nullptr, pkt, outfile);
    /* add sequence end code to have a real MPEG file */
    if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
        fwrite(endcode, 1, sizeof(endcode), outfile);
    fclose(outfile);

    avcodec_free_context(&c);
    av_frame_free(&frame);
    av_packet_free(&pkt);


    sws_freeContext(sws_ctx);
}

void VideoSaver::work() {
    while(true) {
        Optional<Mg::Image2D> job;
        {
            std::unique_lock lock(m_mutex);
            while(m_buffer.empty()){
                if(m_finished) return;
                m_cv.wait(lock);
            }
            job = std::move(m_buffer.front());
            m_buffer.pop();
        }
        encodeSingleFrame(*job);
    }
}

void VideoSaver::appendFrame(Mg::Image2D&& currentFrame) {
    std::lock_guard lock(m_mutex);
    bool wake = m_buffer.empty();
    m_buffer.emplace(std::move(currentFrame));
    if (wake) m_cv.notify_one();
}

void VideoSaver::startRecording(const char* path, const Vector2i& size) {

    const char *codec_name = "libx265";
    int ret;

    /* find the mpeg1video encoder */
    codec = avcodec_find_encoder_by_name(codec_name);
    if (!codec) {
        fprintf(stderr, "Codec '%s' not found\n", codec_name);
        exit(1);
    }

    c = avcodec_alloc_context3(codec);
    if (!c) {
        fprintf(stderr, "Could not allocate video codec context\n");
        exit(1);
    }

    pkt = av_packet_alloc();
    if (!pkt)
        exit(1);

    /* put sample parameters */
    c->bit_rate = 400000;
    /* resolution must be a multiple of two */
    CORRADE_INTERNAL_ASSERT(size.x()%2 == 0 && size.y()%2 == 0);
    c->width = size.x();
    c->height = size.y();
    c->thread_count = 8;

    /* frames per second */
    c->time_base = (AVRational){1, 25};
    c->framerate = (AVRational){25, 1};

    /* emit one intra frame every ten frames
     * check frame pict_type before passing frame
     * to encoder, if frame->pict_type is AV_PICTURE_TYPE_I
     * then gop_size is ignored and the output of encoder
     * will always be I frame irrespective to gop_size
     */
    c->gop_size = 10;
    c->max_b_frames = 1;
    c->pix_fmt = AV_PIX_FMT_YUV420P;

    if (codec->id == AV_CODEC_ID_H264)
        av_opt_set(c->priv_data, "preset", "slow", 0);

    /* open it */
    ret = avcodec_open2(c, codec, nullptr);
    if (ret < 0) {
        Debug{} << "Could not open codec";
        exit(1);
    }

    outfile = fopen(path, "wb");
    if (!outfile) {
        Debug{} << "Could not open" << path;
        exit(1);
    }

    frame = av_frame_alloc();
    if (!frame) {
        Debug{} << "Could not allocate video frame";
        exit(1);
    }
    frame->format = c->pix_fmt;
    frame->width  = c->width;
    frame->height = c->height;

    ret = av_frame_get_buffer(frame, 32);
    if (ret < 0) {
        Debug{} << "Could not allocate the video frame data";
        exit(1);
    }

    /* create scaling context */
    sws_ctx = sws_getContext(size.x(), size.y(), AV_PIX_FMT_RGB24,
                             size.x(), size.y(), AV_PIX_FMT_YUV420P,
                             SWS_BILINEAR, nullptr, nullptr, nullptr);
    if(!sws_ctx) {
        Debug{} << "Could not create scaling context";
    }

    /* allocate source and destination image buffers */
    ret = av_image_alloc(m_data, m_line_size, size.x(), size.y(), AV_PIX_FMT_RGB24, 32);
    if (ret < 0) {
        Debug{} << "Could not allocate source image buffer";
    }

    /* spin up worker thread */
    m_thread = std::thread([this]{ work(); });
}

void VideoSaver::encodeSingleFrame(Mg::Image2D const& image) {
    /* make sure the frame data is writable */
    int ret = av_frame_make_writable(frame);
    if (ret < 0)
        exit(1);

    auto pixel = image.pixels<Color4ub>().transposed<0, 1>().flipped<1>();

    for (size_t y = 0; y < c->height; y++) {
        for (size_t x = 0; x < c->width; x++) {
            m_data[0][y*m_line_size[0] + 3*x] = pixel[x][y].r();
            m_data[0][y*m_line_size[0] + 3*x + 1] = pixel[x][y].g();
            m_data[0][y*m_line_size[0] + 3*x  + 2] = pixel[x][y].b();
        }
    }

    /* convert to destination format */
    sws_scale(sws_ctx, (const uint8_t * const*)m_data,
              m_line_size, 0, frame->height, frame->data, frame->linesize);

    frame->pts = encodedFrames++;

    /* encode the image */
    encode(c, frame, pkt, outfile);
}

}