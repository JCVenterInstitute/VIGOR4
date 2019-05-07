package org.jcvi.vigor.utils;

import org.jcvi.vigor.exception.VigorException;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;

public class Outfiles implements AutoCloseable {

    private final Map<Path, BufferedWriter> writers;
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
            for (BufferedWriter writer : writers.values()) {
                try {
                    writer.close();
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

    public void flush () throws IOException {

        synchronized (writers) {
            List<IOException> exceptions = new ArrayList<>();
            for (BufferedWriter writer : writers.values()) {
                try {
                    writer.flush();
                } catch (IOException e) {
                    exceptions.add(e);
                }
            }
            if (! exceptions.isEmpty()) {
                throw exceptions.get(0);
            }
        }
    }

    /**
     * Get a writer for a given path
     * @param path
     * @return
     */
    public BufferedWriter getWriter(Path path) throws VigorException, IOException {
        Path newPath = getAbsolutePath(path);
        synchronized (writers) {
            if (!writers.containsKey(newPath)) {
                BufferedWriter writer = getBuffer(newPath);
                writers.put(newPath, writer);
            }
            return writers.get(newPath);
        }
    }

    public BufferedWriter getWriter(GenerateVigorOutput.Outfile outfile) throws IOException, VigorException {
        return getWriter(getOutfilePath(outfile));
    }

    private Path getAbsolutePath(Path path) throws VigorException {
        Path newPath =  rootPath.resolve(path).toAbsolutePath();
        if (! newPath.startsWith(rootPath)) {
            throw new VigorException(String.format("%s is not under output directory %s", newPath, rootPath));
        }
        return newPath;
    }

    private Path getOutfilePath(GenerateVigorOutput.Outfile outfile) {
        return rootPath.resolve(baseName + "." + outfile.extension);
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

    private Optional<BufferedWriter> removeWriter(Path path) throws VigorException {
        synchronized (writers) {
            return Optional.ofNullable(writers.get(getAbsolutePath(path)));
        }
    }

    private Optional<BufferedWriter> removeWriter(GenerateVigorOutput.Outfile outfile) throws VigorException {
        return removeWriter(getOutfilePath(outfile));
    }
}
