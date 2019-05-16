package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;
import java.util.function.Consumer;

public class Outfiles implements AutoCloseable {

    class Buffer {

        final BufferedWriter bufferedWriter;
        final Consumer<BufferedWriter> onOpen;
        final Consumer<BufferedWriter> onClose;

        Buffer(BufferedWriter bw, Consumer<BufferedWriter> onOpen, Consumer<BufferedWriter> onClose) {
            bufferedWriter = bw;
            this.onOpen = onOpen;
            this.onClose = onClose;
        }

        void close() throws IOException {
            onClose.accept(bufferedWriter);
            bufferedWriter.close();
        }

        void open() throws IOException {
            onOpen.accept(bufferedWriter);
        }
    }

    private static Logger LOGGER = LogManager.getLogger(Outfiles.class);

    private final Map<Path, Buffer> writers;
    private final Path rootPath;
    private final boolean overwrite;
    private final String baseName;

    public Outfiles(Path rootPath, String baseName, boolean overwrite) {
        this.rootPath = rootPath;
        this.baseName = baseName;
        this.overwrite = overwrite;
        this.writers = new HashMap<>();
    }

    public void close () throws IOException {

        synchronized (writers) {
            List<IOException> exceptions = new ArrayList<>(writers.size());
            for (Map.Entry<Path, Buffer> entry : writers.entrySet()) {
                LOGGER.trace("Closing writer {} for path {}", entry.getValue(), entry.getKey());
                try {
                    entry.getValue().close();
                } catch (IOException e) {
                    exceptions.add(e);
                }
            }
            if (!exceptions.isEmpty()) {
                // TODO join all exceptions somehow meaningfully
                throw exceptions.get(0);
            }
        }
    }

    public void flush() {

        synchronized (writers) {
            List<Path> closed = new ArrayList<>();
            for (Path path: writers.keySet()) {
                LOGGER.trace("Flushing writer for path {}", path);
                try {
                    writers.get(path).bufferedWriter.flush();
                } catch (IOException e) {
                    closed.add(path);
                    LOGGER.debug("flushing path {} got {}:{}", path, e.getClass(), e.getMessage());
                }
            }
            for (Path path: closed) {
                writers.remove(path);
            }
        }
    }

    /**
     * Get a writer for a given path
     * @param path
     * @return
     */
    public BufferedWriter getWriter(Path path, Consumer<BufferedWriter> onOpen, Consumer<BufferedWriter> onClose) throws VigorException, IOException {
        Path newPath = getAbsolutePath(path);
        synchronized (writers) {
            if (!writers.containsKey(newPath)) {
                BufferedWriter writer = getBuffer(newPath);
                Buffer b = new Buffer(writer, onOpen, onClose);
                writers.put(newPath, b);
                b.open();
            }
            return writers.get(newPath).bufferedWriter;
        }
    }

    public BufferedWriter getWriter(Path path) throws VigorException, IOException {
        return getWriter(path, b -> {}, b -> {});
    }

    private Path getAbsolutePath(Path path) throws VigorException {
        Path newPath =  rootPath.resolve(path).toAbsolutePath();
        if (! newPath.startsWith(rootPath)) {
            throw new VigorException(String.format("%s is not under output directory %s", newPath, rootPath));
        }
        return newPath;
    }

    public Path getBaseFilePath(String extension) {
        return rootPath.resolve(baseName + "." + extension);
    }

    private BufferedWriter getBuffer(Path path) throws IOException {
        List<OpenOption> openOptionsList = new ArrayList<>();
        if (overwrite) {
            openOptionsList.add(StandardOpenOption.CREATE);
            openOptionsList.add(StandardOpenOption.TRUNCATE_EXISTING);
        } else {
            openOptionsList.add(StandardOpenOption.CREATE_NEW);
        }
        OpenOption[] openOptions = openOptionsList.toArray(new OpenOption[] {});
        return Files.newBufferedWriter(path, Charset.forName("UTF-8"), openOptions);
    }

    private Optional<Buffer> removeWriter(Path path) throws VigorException {
        synchronized (writers) {
            return Optional.ofNullable(writers.remove(getAbsolutePath(path)));
        }
    }

    public Optional<BufferedWriter> close(Path path) throws IOException, VigorException {
        synchronized (writers) {
            Optional<Buffer> buffer = removeWriter(path);
            if (buffer.isPresent()) {
                Buffer b = buffer.get();
                LOGGER.trace("Closing writer {} for path {}", b.bufferedWriter, path);
                b.close();
            }
            return Optional.ofNullable(buffer == null || ! buffer.isPresent() ? null : buffer.get().bufferedWriter);
        }
    }

}
