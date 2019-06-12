package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class WriterBundle implements AutoCloseable{
    private static Logger LOGGER = LogManager.getLogger(WriterBundle.class);

    private final BufferedWriter[] writers;
    private final Set<BufferedWriter> borrowed = new HashSet<>();

    public WriterBundle(BufferedWriter... writers) {
        this.writers = writers;
    }

    public void write(String value) throws IOException {
        for (BufferedWriter w : writers) {
            try {
                w.write(value);
            } catch (IOException e ) {
                LOGGER.error("problem writing to {}", w);
                throw e;
            }

        }
    }

    public void flush() throws IOException {
        for (BufferedWriter w : writers) {
            w.flush();
        }
    }

    public void write(char c) throws IOException {
        for (BufferedWriter w: writers) {
            try {
                w.write(c);
            } catch (IOException e ) {
                LOGGER.error("problem writing to {}", w);
                throw e;
            }
        }
    }

    public void newLine() throws IOException {
        for (BufferedWriter w: writers) {
            try {
                w.newLine();
            } catch (IOException e ) {
                LOGGER.error("problem writing to {}", w);
                throw e;
            }
        }
    }

    public void borrowed(BufferedWriter writer) {
        LOGGER.trace("{} is borrowed", writer);
        borrowed.add(writer);
    }

    @Override
    public void close() throws IOException {
        LOGGER.trace("borrowed writers: {}", String.join(",", borrowed.stream().map(BufferedWriter::toString).collect(Collectors.toList())));
        for (BufferedWriter writer: writers) {
            if (! borrowed.contains(writer)) {
                LOGGER.trace("Closing {}", writer);
                try {
                    writer.close();
                } catch (IOException e) {
                    LOGGER.warn("error closing writer {} got {}:{}", writer, e.getClass(), e.getMessage());
                }
            }
        }
    }
}
